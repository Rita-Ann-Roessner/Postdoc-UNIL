#!/usr/bin/env python3
"""
Benchmark structure-prediction methods (MD, AF3, ESMFold2, Protenix-v2, OpenDDE)
against ground-truth binder/decoy labels, split by target species (human/mouse).

Produces a 1x5 panel figure: box + strip plot of each method's score,
grouped by species (x-axis) and colored by binder/decoy (hue), annotated
with full AUC and AUC_0.1 (partial AUC restricted to FPR <= 0.1, normalized).

------------------------------------------------------------------------------
EXPECTED INPUT SCHEMA (edit CONFIG below if yours differs)
------------------------------------------------------------------------------
Two CSV files, one per species:
    --human-csv  path/to/human_results.csv
    --mouse-csv  path/to/mouse_results.csv

Each CSV must contain:
    - an identifier column        (default: "name")
    - a ground-truth label column (default: "label", values "binder"/"decoy")
    - one score column per method, e.g.:
        MD:          "MD_dG"        (lower/more negative = stronger binder)
        AF3:         "AF3_ipTM"     (higher = stronger binder)
        ESMFold2:    "ESMFold2_ipTM"
        Protenix-v2: "Protenix-v2_ipTM"
        OpenDDE:     "OpenDDE_ipTM"

If your column names differ, just edit the `METHODS` list in CONFIG below —
no other code changes needed. `higher_is_better=False` flips the AUC/plot
convention for metrics like MD ΔG where more negative = better binder.
------------------------------------------------------------------------------
"""

import argparse
import os
import sys

import inspect

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc as sk_auc

# --- seaborn compatibility -------------------------------------------------
# The `legend=` keyword on categorical plots only exists in seaborn >= 0.13.
# On older seaborn the auto-legend is created per-axes instead, so emulate
# `legend=False` by removing it afterwards.
_SNS_HAS_LEGEND_KW = "legend" in inspect.signature(sns.boxplot).parameters


def _cat_plot(fn, **kwargs):
    want_legend = kwargs.pop("legend", True)
    if _SNS_HAS_LEGEND_KW:
        kwargs["legend"] = want_legend
        return fn(**kwargs)
    # seaborn < 0.13 has no edgecolor="auto" for stripplot
    if kwargs.get("edgecolor") == "auto":
        kwargs["edgecolor"] = "gray"
    ax = fn(**kwargs)
    target = ax if ax is not None else kwargs.get("ax")
    if not want_legend and target is not None and target.get_legend() is not None:
        target.get_legend().remove()
    return ax

# ==============================================================================
# CONFIG — edit this block to match your actual result files
# ==============================================================================

ID_COL = "name"
LABEL_COL = "label"          # values expected: "binder" / "decoy" (case-insensitive)
BINDER_VALUE = "binder"
DECOY_VALUE = "decoy"

# --- Optional: derive LABEL_COL automatically from a BLI-style Classification
# column instead of requiring a pre-computed binder/decoy label per file.
# Used only if --classification-csv is passed on the command line, or if a
# CLASSIFICATION_COL column is already present in the per-species result CSVs.
CLASSIFICATION_ID_COL = "name"   # ID column name inside the classification CSV
CLASSIFICATION_ID_COL_FALLBACKS = ["name", "EPFL Name", "VisterraName"]  # tried in order if the primary is absent
CLASSIFICATION_COL = "Classification"

# species-aware mapping: classification string (lowercased) -> label per species,
# or None to EXCLUDE that row from ground truth for that species (ambiguous /
# no data / control, rather than guessing).
CLASSIFICATION_TO_LABEL = {
    "cross-reactive":                    {"human": "binder", "mouse": "binder"},
    "human-specific":                    {"human": "binder", "mouse": "decoy"},
    "mouse-specific":                    {"human": "decoy",  "mouse": "binder"},
    "non-binding":                       {"human": "decoy",  "mouse": "decoy"},
    "non-binding (no mupd-1 data)":      {"human": "decoy",  "mouse": None},
    "likely nonspecific/sticky":         {"human": "decoy",  "mouse": "decoy"},
    "borderline (weak human trend)":     {"human": None,      "mouse": "decoy"},
    "borderline (weak mouse trend)":     {"human": "decoy",  "mouse": None},
    "positive control (natural ligand)": {"human": None,      "mouse": None},
}

# One entry per panel, in left-to-right plotting order.
METHODS = [
    {"panel_title": "MD",          "column": "MD_dG",             "higher_is_better": False, "ylabel": r"$\Delta G$"},
    {"panel_title": "AF3",         "column": "AF3_ipTM",          "higher_is_better": True,  "ylabel": "ipTM"},
    {"panel_title": "ESMFold2",    "column": "ESMFold2_ipTM",     "higher_is_better": True,  "ylabel": "ipTM"},
    {"panel_title": "Protenix-v2", "column": "Protenix-v2_ipTM",  "higher_is_better": True,  "ylabel": "ipTM"},
    {"panel_title": "OpenDDE",     "column": "OpenDDE_ipTM",      "higher_is_better": True,  "ylabel": "ipTM"},
]

SPECIES_ORDER = ["human", "mouse"]
HUE_ORDER = [BINDER_VALUE, DECOY_VALUE]
PALETTE = {BINDER_VALUE: "#4C72B0", DECOY_VALUE: "#DD8452"}  # matches seaborn 'deep' blue/orange
PARTIAL_FPR_MAX = 0.1  # for AUC_0.1

# --- Cross-reactivity mode: pool human+mouse rows and ask whether a method
# separates genuinely cross-reactive designs from everything else (single-
# species binders + non-binders), rather than splitting by species.
CROSS_REACTIVE_VALUE = "cross-reactive"
SINGLE_NONE_VALUE = "single/none"
CROSS_REACTIVITY_ORDER = [CROSS_REACTIVE_VALUE, SINGLE_NONE_VALUE]
CROSS_REACTIVITY_PALETTE = {CROSS_REACTIVE_VALUE: "#4C72B0", SINGLE_NONE_VALUE: "#DD8452"}
CROSS_REACTIVITY_TICK_LABELS = ["cross-\nreactive", "single/\nnone"]

# Top-fraction cutoffs for Enrichment Factor / Precision-at-K (prevalence-aware
# complements to AUC, which is prevalence-invariant by construction).
ENRICHMENT_TOP_FRACS = [0.10, 0.20]

# ==============================================================================


def compute_auc_metrics(y_true_binder, scores, higher_is_better):
    """
    y_true_binder: array-like of bool/0-1, True/1 = binder (positive class)
    scores: array-like of the raw metric values
    higher_is_better: if False, sign-flip scores so that "higher = more binder-like"
                       before computing ROC (needed for e.g. MD deltaG).
    Returns (auc_full, auc_partial_normalized) or (nan, nan) if degenerate
    (e.g. only one class present, or too few points).
    """
    y_true_binder = np.asarray(y_true_binder).astype(int)
    scores = np.asarray(scores, dtype=float)

    mask = ~np.isnan(scores)
    y_true_binder = y_true_binder[mask]
    scores = scores[mask]

    if len(np.unique(y_true_binder)) < 2 or len(scores) < 2:
        return np.nan, np.nan

    score_for_roc = scores if higher_is_better else -scores

    fpr, tpr, _ = roc_curve(y_true_binder, score_for_roc)
    auc_full = sk_auc(fpr, tpr)

    # Partial AUC restricted to FPR in [0, PARTIAL_FPR_MAX], normalized by
    # PARTIAL_FPR_MAX so the metric is comparable to full AUC (max = 1.0).
    # Interpolate tpr at exactly PARTIAL_FPR_MAX for a clean boundary.
    if PARTIAL_FPR_MAX >= fpr.max():
        auc_partial_raw = auc_full * min(PARTIAL_FPR_MAX, fpr.max())
    else:
        # points strictly within the window
        within = fpr <= PARTIAL_FPR_MAX
        fpr_w = fpr[within].tolist()
        tpr_w = tpr[within].tolist()
        # interpolate the boundary point at PARTIAL_FPR_MAX
        tpr_at_cutoff = np.interp(PARTIAL_FPR_MAX, fpr, tpr)
        fpr_w.append(PARTIAL_FPR_MAX)
        tpr_w.append(tpr_at_cutoff)
        trapz_fn = getattr(np, "trapezoid", None) or np.trapz
        auc_partial_raw = trapz_fn(tpr_w, fpr_w)

    auc_partial_norm = auc_partial_raw / PARTIAL_FPR_MAX
    return auc_full, auc_partial_norm


def compute_enrichment_and_precision(y_true_pos, scores, higher_is_better, top_fracs=ENRICHMENT_TOP_FRACS):
    """
    Prevalence-aware complements to AUC. AUC/AUC_0.1 are prevalence-invariant
    by construction (they depend only on the score distribution within each
    class, not on how many of each class there are) -- useful for comparing
    across datasets, but they don't tell you how practically useful the
    ranking is when positives are rare.

    For each top_frac X (e.g. 0.10 = top 10% of the ranked set by score):
      - precision_topX  = fraction of the top-X% that are true positives
      - EF_topX          = precision_topX / prevalence
                            (how many times better than random selection)
    prevalence = overall fraction of positives in the full set.

    Returns a dict, e.g. {"prevalence": 0.10, "precision_top10": 0.4,
    "EF_top10": 4.0, "precision_top20": 0.3, "EF_top20": 3.0}. NaNs if
    degenerate (no positives, or too few points).
    """
    y_true_pos = np.asarray(y_true_pos).astype(int)
    scores = np.asarray(scores, dtype=float)

    mask = ~np.isnan(scores)
    y_true_pos = y_true_pos[mask]
    scores = scores[mask]

    n = len(y_true_pos)
    result = {}
    if n == 0 or y_true_pos.sum() == 0:
        result["prevalence"] = np.nan if n == 0 else 0.0
        for frac in top_fracs:
            pct = int(round(frac * 100))
            result[f"precision_top{pct}"] = np.nan
            result[f"EF_top{pct}"] = np.nan
        return result

    prevalence = y_true_pos.sum() / n
    result["prevalence"] = prevalence

    score_for_ranking = scores if higher_is_better else -scores
    order = np.argsort(-score_for_ranking)  # descending: best first
    sorted_labels = y_true_pos[order]

    for frac in top_fracs:
        pct = int(round(frac * 100))
        k = max(1, int(round(frac * n)))
        precision_k = sorted_labels[:k].mean()
        ef_k = precision_k / prevalence if prevalence > 0 else np.nan
        result[f"precision_top{pct}"] = precision_k
        result[f"EF_top{pct}"] = ef_k

    return result


def derive_label_from_classification(classification, species):
    """Map a BLI Classification string to 'binder'/'decoy'/None for a given species."""
    key = str(classification).strip().lower()
    entry = CLASSIFICATION_TO_LABEL.get(key)
    if entry is None:
        return None  # unrecognized category -> exclude rather than guess
    return entry.get(species)


def load_species_df(csv_path, species_name, classification_df=None):
    df = pd.read_csv(csv_path)
    if ID_COL not in df.columns:
        raise ValueError(
            f"{csv_path}: missing required ID column '{ID_COL}'. "
            f"Found columns: {list(df.columns)}"
        )
    df = df.copy()

    # Merge in an external classification table if provided. If the score CSV
    # already carries the classification column (the merged tables written by
    # Analysis/build_dataset.py do), skip the merge -- otherwise pandas would
    # silently rename both sides to Classification_x / Classification_y.
    if classification_df is not None and CLASSIFICATION_COL not in df.columns:
        df = df.merge(classification_df, on=ID_COL, how="left")
        n_unmatched = df[CLASSIFICATION_COL].isna().sum()
        if n_unmatched:
            print(
                f"  [warn] {csv_path}: {n_unmatched} row(s) had no matching "
                f"'{ID_COL}' in the classification CSV and will be excluded."
            )

    if LABEL_COL not in df.columns:
        if CLASSIFICATION_COL not in df.columns:
            raise ValueError(
                f"{csv_path}: no '{LABEL_COL}' column and no '{CLASSIFICATION_COL}' "
                f"column (directly or via --classification-csv) to derive it from."
            )
        df[LABEL_COL] = df[CLASSIFICATION_COL].apply(
            lambda c: derive_label_from_classification(c, species_name)
        )
        n_before = len(df)
        excluded_counts = (
            df.loc[df[LABEL_COL].isna(), CLASSIFICATION_COL]
            .value_counts()
            .to_dict()
        )
        df = df.dropna(subset=[LABEL_COL])
        n_excluded = n_before - len(df)
        if n_excluded:
            print(
                f"  [{species_name}] excluded {n_excluded}/{n_before} row(s) as "
                f"ambiguous/no-data/control for this species: {excluded_counts}"
            )
    else:
        df[LABEL_COL] = df[LABEL_COL].astype(str).str.strip().str.lower()
        bad_labels = set(df[LABEL_COL].unique()) - {BINDER_VALUE, DECOY_VALUE}
        if bad_labels:
            raise ValueError(
                f"{csv_path}: unexpected label value(s) {bad_labels}; "
                f"expected only '{BINDER_VALUE}'/'{DECOY_VALUE}'."
            )

    # The natural ligand (HuPDL1-Fc) is labelled "Cross-reactive" in the BLI
    # sheet but is a positive control, not a de novo design: 227 aa and a
    # PDB-deposited complex. Drop it by name so it cannot inflate the AUCs.
    n_before = len(df)
    df = df[df[ID_COL] != "HuPDL1-Fc"]
    if len(df) < n_before:
        print(f"  [{species_name}] excluded the HuPDL1-Fc positive control "
              f"(natural ligand, not a design).")

    df["PDL1"] = species_name
    return df


def _derive_out_paths(base_out, base_out_csv, mode_suffix):
    """Insert a mode suffix before the file extension, e.g.
    'method_benchmark.png' -> 'method_benchmark_species.png'."""
    root, ext = os.path.splitext(base_out)
    out_path = f"{root}_{mode_suffix}{ext}"
    out_csv_path = None
    if base_out_csv:
        croot, cext = os.path.splitext(base_out_csv)
        out_csv_path = f"{croot}_{mode_suffix}{cext}"
    return out_path, out_csv_path


def build_per_design_cross_reactivity_table(human_df, mouse_df, classification_df, methods_active):
    """
    Collapse human+mouse rows into ONE row per design per method: the
    "worst-case" (bottleneck) score across the two species, since a design
    can only be truly cross-reactive if it performs adequately against BOTH
    targets -- its weaker side is what would sink it in a real selection.

    For higher_is_better methods (ipTM): combined_score = min(human, mouse).
    For lower_is_better methods (MD dG): combined_score = max(human, mouse)
    (max ΔG = the less favorable/weaker of the two).

    Returns a long-format DataFrame: [name, method, combined_score, CrossReactivity]
    Designs missing a valid score in EITHER species (for a given method) are
    dropped for that method, with a printed count.
    """
    class_lookup = classification_df.drop_duplicates(subset=[ID_COL]).copy()
    class_lookup[CLASSIFICATION_COL] = class_lookup[CLASSIFICATION_COL].astype(str).str.strip()

    records = []
    for m in methods_active:
        col = m["column"]
        carry = [ID_COL] + ([CLASSIFICATION_COL] if CLASSIFICATION_COL in human_df.columns else [])
        h = human_df[carry + [col]].rename(columns={col: "human_score"})
        mo = mouse_df[[ID_COL, col]].rename(columns={col: "mouse_score"})
        merged = h.merge(mo, on=ID_COL, how="inner")

        n_before = len(merged)
        merged = merged.dropna(subset=["human_score", "mouse_score"])
        n_after = len(merged)
        if n_after < n_before:
            print(f"  [{m['panel_title']}] {n_before - n_after} design(s) missing a score in "
                  f"one species -- excluded from per-design cross-reactivity analysis.")

        if m["higher_is_better"]:
            merged["combined_score"] = np.minimum(merged["human_score"], merged["mouse_score"])
        else:
            merged["combined_score"] = np.maximum(merged["human_score"], merged["mouse_score"])

        if CLASSIFICATION_COL not in merged.columns:
            merged = merged.merge(class_lookup[[ID_COL, CLASSIFICATION_COL]],
                                  on=ID_COL, how="left")
        n_unmatched = merged[CLASSIFICATION_COL].isna().sum()
        if n_unmatched:
            print(f"  [{m['panel_title']}] {n_unmatched} design(s) had no matching classification "
                  f"and will be excluded.")

        is_positive_control = (
            (merged[CLASSIFICATION_COL].str.lower() == "positive control (natural ligand)")
            | (merged[ID_COL] == "HuPDL1-Fc")
        )
        merged = merged[~is_positive_control]

        is_cross_reactive = merged[CLASSIFICATION_COL].str.lower() == "cross-reactive"
        merged["CrossReactivity"] = np.where(is_cross_reactive, CROSS_REACTIVE_VALUE, SINGLE_NONE_VALUE)
        merged.loc[merged[CLASSIFICATION_COL].isna(), "CrossReactivity"] = np.nan
        merged = merged.dropna(subset=["CrossReactivity"])

        merged["method"] = m["panel_title"]
        records.append(merged[[ID_COL, "method", "combined_score", "CrossReactivity"]])

    return pd.concat(records, ignore_index=True) if records else pd.DataFrame(
        columns=[ID_COL, "method", "combined_score", "CrossReactivity"]
    )


def render_species_mode(combined, methods_active, out_path, out_csv_path, dpi):
    n_panels = len(methods_active)
    fig, axes = plt.subplots(1, n_panels, figsize=(3.2 * n_panels, 4.2), sharey=False)
    if n_panels == 1:
        axes = [axes]

    auc_records = []

    for ax, m in zip(axes, methods_active):
        col = m["column"]
        panel_df = combined[[ID_COL, "PDL1", LABEL_COL, col]].dropna(subset=[col])

        _cat_plot(
            sns.boxplot,
            data=panel_df, x="PDL1", y=col, hue=LABEL_COL,
            order=SPECIES_ORDER, hue_order=HUE_ORDER, palette=PALETTE,
            ax=ax, showfliers=False, width=0.6,
            boxprops=dict(alpha=0.85), legend=(ax is axes[-1]),
        )
        _cat_plot(
            sns.stripplot,
            data=panel_df, x="PDL1", y=col, hue=LABEL_COL,
            order=SPECIES_ORDER, hue_order=HUE_ORDER, palette=PALETTE,
            ax=ax, dodge=True, size=4, alpha=0.6, linewidth=0.4,
            edgecolor="auto", legend=False,
        )

        # AUC / AUC_0.1 per species, annotated at the top of each species column
        y_top = panel_df[col].max()
        y_range = panel_df[col].max() - panel_df[col].min()
        y_range = y_range if y_range > 0 else 1.0

        for xi, species in enumerate(SPECIES_ORDER):
            sub = panel_df[panel_df["PDL1"] == species]
            y_true_binder = (sub[LABEL_COL] == BINDER_VALUE).values
            auc_full, auc_partial = compute_auc_metrics(
                y_true_binder, sub[col].values, m["higher_is_better"]
            )
            enrich = compute_enrichment_and_precision(
                y_true_binder, sub[col].values, m["higher_is_better"]
            )
            auc_records.append({
                "method": m["panel_title"], "species": species,
                "AUC": auc_full, "AUC_0.1": auc_partial,
                "n_binder": int(y_true_binder.sum()), "n_decoy": int((~y_true_binder).sum()),
                **enrich,
            })
            label_txt = f"AUC={auc_full:.2f}\n" + r"AUC$_{0.1}$" + f"={auc_partial:.2f}"
            ax.text(
                xi, y_top + 0.06 * y_range, label_txt,
                ha="center", va="bottom", fontsize=8.5,
            )

        ax.set_title(m["panel_title"])
        ax.set_xlabel("PDL1" if ax is axes[len(axes) // 2] else "")
        ax.set_ylabel(m["ylabel"])
        # headroom for the annotation text
        ax.set_ylim(panel_df[col].min() - 0.08 * y_range, y_top + 0.28 * y_range)
        if m["panel_title"] == "MD":
            ax.invert_yaxis()  # more negative (better) drawn upward, matching example figure

    # single shared legend on the last panel only
    handles, labels = axes[-1].get_legend_handles_labels()
    if handles:
        axes[-1].legend(handles[:2], labels[:2], loc="lower right", frameon=True, fontsize=8, title=None)

    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    print(f"Wrote figure: {out_path}")

    auc_df = pd.DataFrame(auc_records)
    print(f"\n=== AUC SUMMARY (species) ===")
    print(auc_df.to_string(index=False))

    if out_csv_path:
        auc_df.to_csv(out_csv_path, index=False)
        print(f"Wrote AUC table: {out_csv_path}")
    plt.close(fig)


def render_cross_reactivity_mode(human_df, mouse_df, classification_df, methods_active,
                                  out_path, out_csv_path, dpi):
    """
    One box per method: 'cross-reactive' vs 'single/none', using ONE combined
    (worst-case, bottleneck) score per design -- this reflects the actual task
    of picking cross-reactive candidates out of the full set, rather than
    treating each design's human and mouse scores as independent samples.
    """
    long_df = build_per_design_cross_reactivity_table(human_df, mouse_df, classification_df, methods_active)

    n_panels = len(methods_active)
    fig, axes = plt.subplots(1, n_panels, figsize=(3.2 * n_panels, 4.2), sharey=False)
    if n_panels == 1:
        axes = [axes]

    auc_records = []

    for ax, m in zip(axes, methods_active):
        panel_df = long_df[long_df["method"] == m["panel_title"]]
        if panel_df.empty:
            ax.set_title(m["panel_title"])
            ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center", va="center")
            continue

        _cat_plot(
            sns.boxplot,
            data=panel_df, x="CrossReactivity", y="combined_score", hue="CrossReactivity",
            order=CROSS_REACTIVITY_ORDER, hue_order=CROSS_REACTIVITY_ORDER,
            palette=CROSS_REACTIVITY_PALETTE,
            ax=ax, showfliers=False, width=0.6,
            boxprops=dict(alpha=0.85), legend=False, dodge=False,
        )
        _cat_plot(
            sns.stripplot,
            data=panel_df, x="CrossReactivity", y="combined_score", hue="CrossReactivity",
            order=CROSS_REACTIVITY_ORDER, hue_order=CROSS_REACTIVITY_ORDER,
            palette=CROSS_REACTIVITY_PALETTE,
            ax=ax, dodge=False, size=4, alpha=0.6, linewidth=0.4,
            edgecolor="auto", legend=False,
        )

        y_top = panel_df["combined_score"].max()
        y_range = panel_df["combined_score"].max() - panel_df["combined_score"].min()
        y_range = y_range if y_range > 0 else 1.0

        y_true_cross = (panel_df["CrossReactivity"] == CROSS_REACTIVE_VALUE).values
        auc_full, auc_partial = compute_auc_metrics(
            y_true_cross, panel_df["combined_score"].values, m["higher_is_better"]
        )
        enrich = compute_enrichment_and_precision(
            y_true_cross, panel_df["combined_score"].values, m["higher_is_better"]
        )
        auc_records.append({
            "method": m["panel_title"],
            "AUC": auc_full, "AUC_0.1": auc_partial,
            "n_cross_reactive": int(y_true_cross.sum()),
            "n_single_none": int((~y_true_cross).sum()),
            **enrich,
        })
        label_txt = f"AUC={auc_full:.2f}\n" + r"AUC$_{0.1}$" + f"={auc_partial:.2f}"
        ax.text(
            0.5, 0.98, label_txt, transform=ax.transAxes,
            ha="center", va="top", fontsize=8.5,
        )

        ax.set_title(m["panel_title"])
        ax.set_xlabel("PDL1" if ax is axes[len(axes) // 2] else "")
        ax.set_ylabel(m["ylabel"])
        ax.set_xticks([0, 1])
        ax.set_xticklabels(CROSS_REACTIVITY_TICK_LABELS)
        ax.set_ylim(panel_df["combined_score"].min() - 0.08 * y_range, y_top + 0.20 * y_range)
        if m["panel_title"] == "MD":
            ax.invert_yaxis()

    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    print(f"Wrote figure: {out_path}")

    auc_df = pd.DataFrame(auc_records)
    print(f"\n=== AUC SUMMARY (cross_reactivity, per-design worst-case score) ===")
    print(auc_df.to_string(index=False))

    if out_csv_path:
        auc_df.to_csv(out_csv_path, index=False)
        print(f"Wrote AUC table: {out_csv_path}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--human-csv", required=True, help="Path to human results CSV")
    parser.add_argument("--mouse-csv", required=True, help="Path to mouse results CSV")
    parser.add_argument(
        "--classification-csv", default=None,
        help=(
            "Optional path to a BLI-style CSV with an ID column "
            f"('{CLASSIFICATION_ID_COL}') and a '{CLASSIFICATION_COL}' column "
            "(e.g. 'Cross-reactive', 'Human-specific', ...). If given, binder/decoy "
            "labels are derived per-species automatically (see CLASSIFICATION_TO_LABEL "
            "in the CONFIG block) instead of requiring a pre-computed 'label' column."
        ),
    )
    parser.add_argument(
        "--out", default="method_benchmark.png",
        help="Output figure path (png/pdf/svg). With --mode both, a suffix is "
             "inserted before the extension for each figure (e.g. "
             "'method_benchmark_species.png', 'method_benchmark_cross_reactivity.png').",
    )
    parser.add_argument(
        "--out-csv", default=None,
        help="Optional: base path for the AUC/AUC_0.1 table CSV(s); same suffixing as --out with --mode both.",
    )
    parser.add_argument("--dpi", type=int, default=150)
    parser.add_argument(
        "--mode", choices=["species", "cross_reactivity", "both"], default="species",
        help=(
            "'species' (default): split panels by human/mouse, hue by binder/decoy. "
            "'cross_reactivity': pool human+mouse rows, single box per method for "
            "'cross-reactive' vs 'single/none' (requires --classification-csv). "
            "'both': generate both figures/tables in one run (requires --classification-csv)."
        ),
    )
    args = parser.parse_args()

    modes_to_run = ["species", "cross_reactivity"] if args.mode == "both" else [args.mode]

    if any(m == "cross_reactivity" for m in modes_to_run) and not args.classification_csv:
        raise ValueError("--mode cross_reactivity (or 'both') requires --classification-csv to be provided.")

    classification_df = None
    if args.classification_csv:
        classification_df = pd.read_csv(args.classification_csv)

        resolved_id_col = None
        for candidate in [CLASSIFICATION_ID_COL] + CLASSIFICATION_ID_COL_FALLBACKS:
            if candidate in classification_df.columns:
                resolved_id_col = candidate
                break
        if resolved_id_col is None:
            raise ValueError(
                f"{args.classification_csv}: none of the expected ID columns "
                f"{[CLASSIFICATION_ID_COL] + CLASSIFICATION_ID_COL_FALLBACKS} found. "
                f"Found columns: {list(classification_df.columns)}"
            )
        if resolved_id_col != CLASSIFICATION_ID_COL:
            print(f"  [info] using '{resolved_id_col}' as the classification-CSV ID column "
                  f"(primary '{CLASSIFICATION_ID_COL}' not found).")

        if CLASSIFICATION_COL not in classification_df.columns:
            raise ValueError(
                f"{args.classification_csv}: missing '{CLASSIFICATION_COL}' column. "
                f"Found columns: {list(classification_df.columns)}"
            )
        classification_df = classification_df[[resolved_id_col, CLASSIFICATION_COL]].rename(
            columns={resolved_id_col: ID_COL}
        )

    human_df = load_species_df(args.human_csv, "human", classification_df)
    mouse_df = load_species_df(args.mouse_csv, "mouse", classification_df)

    # Only plot methods whose score column is actually present in BOTH species'
    # CSVs. This lets the script run as more methods get added over time,
    # instead of erroring out or requiring you to hand-edit METHODS.
    available_methods = []
    skipped_methods = []
    for m in METHODS:
        if m["column"] in human_df.columns and m["column"] in mouse_df.columns:
            available_methods.append(m)
        else:
            skipped_methods.append(m["panel_title"])
    if skipped_methods:
        print(f"Skipping method(s) not yet present in both CSVs: {skipped_methods}")
    if not available_methods:
        raise ValueError(
            "None of the configured METHODS columns were found in the provided CSVs. "
            f"Configured columns: {[m['column'] for m in METHODS]}. "
            f"human_df columns: {list(human_df.columns)}. "
            f"mouse_df columns: {list(mouse_df.columns)}."
        )
    METHODS_ACTIVE = available_methods

    combined = pd.concat([human_df, mouse_df], ignore_index=True)

    # Restrict to complete cases: only designs with a valid score for every
    # *available* method, so every panel/AUC comparison uses the exact same
    # set of candidates.
    method_cols = [m["column"] for m in METHODS_ACTIVE]
    n_before = len(combined)
    combined = combined.dropna(subset=method_cols)
    n_after = len(combined)
    if n_after < n_before:
        print(f"Restricting to complete cases (all {len(method_cols)} available method score(s) present): "
              f"{n_after}/{n_before} rows kept.")
    if n_after == 0:
        raise ValueError(
            f"No rows have a valid score for all {len(method_cols)} available method(s) "
            "simultaneously — nothing to plot. Check for column-name mismatches or sparse coverage."
        )

    for mode in modes_to_run:
        if args.mode == "both":
            out_path, out_csv_path = _derive_out_paths(args.out, args.out_csv, mode)
        else:
            out_path, out_csv_path = args.out, args.out_csv

        if mode == "species":
            render_species_mode(combined, METHODS_ACTIVE, out_path, out_csv_path, args.dpi)
        else:
            # Use FRESH, unfiltered reads of the score CSVs here (not the
            # `combined`/`human_df`/`mouse_df` used for species mode, which
            # already dropped borderline/no-data rows per-species). The
            # per-design cross-reactivity table needs both raw scores so it
            # can make its own inclusion/exclusion decision at the design level.
            human_raw = pd.read_csv(args.human_csv)
            mouse_raw = pd.read_csv(args.mouse_csv)
            render_cross_reactivity_mode(
                human_raw, mouse_raw, classification_df, METHODS_ACTIVE,
                out_path, out_csv_path, args.dpi,
            )


if __name__ == "__main__":
    main()