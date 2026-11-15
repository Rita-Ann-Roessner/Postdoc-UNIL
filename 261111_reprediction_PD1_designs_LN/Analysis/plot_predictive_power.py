#!/usr/bin/env python3
"""
Figures for the ESMFold2 / Protenix-v2 predictive-power benchmark.

Run Analysis/predictive_power.py first (or just run this -- it re-runs the
analysis in-process and reuses the same score vectors, so figures and tables
can never drift apart).

  fig1_score_distributions.png  raw score distributions, binder vs decoy
  fig2_roc.png                  ROC curves, each method alone and fused
  fig3_auc_forest.png           AUC +/- 95% CI for every predictor, all 3 tasks
  fig4_method_agreement.png     ESMFold2 vs Protenix scores, quadrant analysis
  fig5_screening_yield.png      what you actually get if you order the top X%
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy import stats
from sklearn.metrics import roc_curve, roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import predictive_power as pp   # noqa: E402

FIG = pp.FIG
E, P = pp.E, pp.P

TASK_LABEL = {"human": "hPD-1 binding", "mouse": "mPD-1 binding",
              "cross_reactive": "cross-reactivity\n(worst-of-both-species score)"}
TASK_ORDER = ["human", "mouse", "cross_reactive"]

C_BINDER, C_DECOY = "#4C72B0", "#DD8452"
GROUP_COLOR = {
    f"{E} alone":          "#4C72B0",
    f"{P} alone":          "#DD8452",
    "Together (fusion)":   "#55A868",
    "Together (learned)":  "#2E7A4F",
    "Together (consensus)": "#8172B2",
    "Baseline":            "#999999",
}

# curated, fixed row order so the three task panels stay row-comparable
FOREST_ROWS = [
    (f"{E} interface pLDDT", f"{E} alone"),
    (f"{E} ipTM", f"{E} alone"),
    (f"{E} pTM", f"{E} alone"),
    (f"{E} pLDDT", f"{E} alone"),
    (f"{E} pAE", f"{E} alone"),
    (f"{E} composite (CV)", f"{E} alone"),
    (f"{P} interface pLDDT", f"{P} alone"),
    (f"{P} ipTM", f"{P} alone"),
    (f"{P} pTM", f"{P} alone"),
    (f"{P} pLDDT", f"{P} alone"),
    (f"{P} i-pAE", f"{P} alone"),
    (f"{P} composite (CV)", f"{P} alone"),
    ("mean-z(ipTM)", "Together (fusion)"),
    ("min-z(ipTM)", "Together (fusion)"),
    ("max-z(ipTM)", "Together (fusion)"),
    ("rank-mean(ipTM)", "Together (fusion)"),
    ("Both ipTM (CV)", "Together (learned)"),
    ("Both ipTM + consensus (CV)", "Together (learned)"),
    ("All metrics, both methods (CV)", "Together (learned)"),
    ("cross-method pose RMSD", "Together (consensus)"),
    ("cross-method epitope overlap", "Together (consensus)"),
    ("binder length (baseline)", "Baseline"),
]

ROC_CURVES = [
    (f"{E} ipTM", "#4C72B0", "-"),
    (f"{P} ipTM", "#DD8452", "-"),
    ("mean-z(ipTM)", "#55A868", "-"),
    ("Both ipTM (CV)", "#2E7A4F", "--"),
    (f"{E} interface pLDDT", "#4C72B0", ":"),
]


def get_scores(store, task):
    df, preds, cvp = store[task]
    allp = {k: v[0] for k, v in preds.items()}
    allp.update({k: v[0] for k, v in cvp.items()})
    return df, allp


# ------------------------------------------------------------------ figure 1
def fig1(store):
    panels = [(f"{E} ipTM", "ipTM"), (f"{P} ipTM", "ipTM"),
              ("mean-z(ipTM)", "fused z-score"),
              (f"{E} interface pLDDT", "interface pLDDT")]
    fig, axes = plt.subplots(len(TASK_ORDER), len(panels),
                             figsize=(3.1 * len(panels), 3.0 * len(TASK_ORDER)))
    for i, task in enumerate(TASK_ORDER):
        df, allp = get_scores(store, task)
        y = df["label"].values.astype(int)
        for j, (key, ylab) in enumerate(panels):
            ax = axes[i, j]
            s = allp[key]
            groups = [s[y == 1], s[y == 0]]
            bp = ax.boxplot(groups, positions=[0, 1], widths=0.55,
                            showfliers=False, patch_artist=True,
                            medianprops=dict(color="k", lw=1.3))
            for patch, c in zip(bp["boxes"], [C_BINDER, C_DECOY]):
                patch.set_facecolor(c)
                patch.set_alpha(0.55)
            rng = np.random.default_rng(1)
            for xi, (vals, c) in enumerate(zip(groups, [C_BINDER, C_DECOY])):
                ax.scatter(xi + rng.uniform(-0.16, 0.16, len(vals)), vals,
                           s=16, color=c, alpha=0.75, linewidth=0.4,
                           edgecolor="white", zorder=3)
            a = roc_auc_score(y, s)
            lo, hi = pp.bootstrap_auc_ci(y, s, n_boot=2000)
            u = stats.mannwhitneyu(groups[0], groups[1], alternative="greater")[1]
            ax.set_title(f"AUC = {a:.2f} [{lo:.2f}-{hi:.2f}]   p = {u:.1g}",
                         fontsize=9)
            ax.set_xticks([0, 1])
            ax.set_xticklabels([f"binder\n(n={int(y.sum())})",
                                f"non-binder\n(n={int((1-y).sum())})"], fontsize=8)
            if j == 0:
                ax.set_ylabel(f"{TASK_LABEL[task]}\n\n{ylab}", fontsize=9)
            else:
                ax.set_ylabel(ylab, fontsize=9)
            if i == 0:
                ax.text(0.5, 1.22, key, transform=ax.transAxes, ha="center",
                        fontsize=10.5, fontweight="bold")
            ax.grid(axis="y", alpha=0.25, lw=0.6)
            ax.set_axisbelow(True)
    fig.suptitle("Score distributions by BLI outcome  (95% CI from stratified "
                 "bootstrap; p from one-sided Mann-Whitney U)",
                 y=1.005, fontsize=11)
    fig.tight_layout()
    out = os.path.join(FIG, "fig1_score_distributions.png")
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("wrote", os.path.relpath(out, pp.ROOT))


# ------------------------------------------------------------------ figure 2
def fig2(store):
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.5))
    for ax, task in zip(axes, TASK_ORDER):
        df, allp = get_scores(store, task)
        y = df["label"].values.astype(int)
        for key, color, ls in ROC_CURVES:
            if key not in allp:
                continue
            s = allp[key]
            fpr, tpr, _ = roc_curve(y, s)
            ax.plot(fpr, tpr, color=color, ls=ls, lw=1.9,
                    label=f"{key}  ({roc_auc_score(y, s):.2f})")
        ax.plot([0, 1], [0, 1], color="k", lw=0.9, ls="--", alpha=0.5,
                label="random  (0.50)")
        ax.axvspan(0, pp.PARTIAL_FPR_MAX, color="grey", alpha=0.10, lw=0)
        ax.text(pp.PARTIAL_FPR_MAX / 2, 0.02, "AUC$_{0.1}$\nregion",
                fontsize=7, ha="center", va="bottom", color="grey")
        ax.set_title(f"{TASK_LABEL[task]}\n"
                     f"{int(y.sum())} positives / {len(y)} designs", fontsize=10)
        ax.set_xlabel("false-positive rate")
        ax.set_ylabel("true-positive rate")
        ax.set_xlim(-0.01, 1.01)
        ax.set_ylim(-0.01, 1.01)
        ax.legend(loc="lower right", fontsize=7.5, frameon=True)
        ax.grid(alpha=0.25, lw=0.6)
        ax.set_axisbelow(True)
    fig.suptitle("ROC: each method alone vs the two together", y=1.02, fontsize=12)
    fig.tight_layout()
    out = os.path.join(FIG, "fig2_roc.png")
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("wrote", os.path.relpath(out, pp.ROOT))


# ------------------------------------------------------------------ figure 3
def fig3(auc_df):
    rows = [r for r in FOREST_ROWS
            if ((auc_df.predictor == r[0]).any())]
    fig, axes = plt.subplots(1, 3, figsize=(15, 0.34 * len(rows) + 2.4),
                             sharey=True)
    ypos = np.arange(len(rows))[::-1]
    for ax, task in zip(axes, TASK_ORDER):
        sub = auc_df[auc_df.task == task].set_index("predictor")
        for yv, (name, group) in zip(ypos, rows):
            if name not in sub.index:
                continue
            r = sub.loc[name]
            c = GROUP_COLOR[group]
            ax.plot([r["AUC_lo"], r["AUC_hi"]], [yv, yv], color=c, lw=2.0,
                    alpha=0.75, solid_capstyle="round")
            ax.scatter([r["AUC"]], [yv], color=c, s=46, zorder=4,
                       marker="D" if r["estimate"] == "out-of-fold" else "o",
                       edgecolor="white", linewidth=0.7)
        ax.axvline(0.5, color="k", ls="--", lw=0.9, alpha=0.6)
        best = sub.loc[[n for n, _ in rows if n in sub.index], "AUC"].max()
        ax.axvline(best, color="#333333", ls=":", lw=0.9, alpha=0.6)
        ax.set_title(f"{TASK_LABEL[task]}\n"
                     f"{int(sub['n_pos'].iloc[0])} positives / "
                     f"{int(sub['n'].iloc[0])} designs", fontsize=10)
        ax.set_xlabel("AUC (95% CI)")
        ax.set_xlim(0.15, 1.0)
        ax.grid(axis="x", alpha=0.25, lw=0.6)
        ax.set_axisbelow(True)
    axes[0].set_yticks(ypos)
    axes[0].set_yticklabels([n for n, _ in rows], fontsize=8.5)
    handles = [Patch(facecolor=c, label=g) for g, c in GROUP_COLOR.items()]
    handles += [Line2D([], [], marker="o", ls="", color="k", label="apparent AUC"),
                Line2D([], [], marker="D", ls="", color="k",
                       label="out-of-fold AUC (20x5-fold CV)")]
    axes[-1].legend(handles=handles, loc="lower right", fontsize=7.5,
                    frameon=True)
    fig.suptitle("Predictive power of every candidate score, alone and together",
                 y=1.0, fontsize=12.5)
    fig.tight_layout()
    out = os.path.join(FIG, "fig3_auc_forest.png")
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("wrote", os.path.relpath(out, pp.ROOT))


# ------------------------------------------------------------------ figure 4
def fig4(wide, tasks):
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.8))
    for ax, task in zip(axes, TASK_ORDER):
        df = tasks[task]
        y = df["label"].values.astype(int)
        x, z = df[f"{E}_ipTM"].values, df[f"{P}_ipTM"].values
        rho = stats.spearmanr(x, z)[0]
        mx, mz = np.median(x), np.median(z)
        ax.axvline(mx, color="grey", lw=0.9, ls="--", alpha=0.7)
        ax.axhline(mz, color="grey", lw=0.9, ls="--", alpha=0.7)
        ax.scatter(x[y == 0], z[y == 0], s=26, c=C_DECOY, alpha=0.7,
                   edgecolor="white", linewidth=0.5, label="non-binder")
        ax.scatter(x[y == 1], z[y == 1], s=54, c=C_BINDER, alpha=0.95,
                   edgecolor="k", linewidth=0.6, label="binder", zorder=3)
        # quadrant precision
        for qx, qz, ha, va in [(1, 1, "right", "top"), (0, 1, "left", "top"),
                               (1, 0, "right", "bottom"), (0, 0, "left", "bottom")]:
            m = ((x >= mx) if qx else (x < mx)) & ((z >= mz) if qz else (z < mz))
            if m.sum() == 0:
                continue
            ax.text(0.97 if qx else 0.03, 0.97 if qz else 0.03,
                    f"{int(y[m].sum())}/{int(m.sum())}\n{y[m].mean():.0%} binders",
                    transform=ax.transAxes, ha=ha, va=va, fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="grey",
                              alpha=0.85))
        ax.set_xlabel(f"{E} ipTM" + ("  (worst species)" if task == "cross_reactive" else ""))
        ax.set_ylabel(f"{P} ipTM" + ("  (worst species)" if task == "cross_reactive" else ""))
        ax.set_title(f"{TASK_LABEL[task]}\nSpearman $\\rho$ = {rho:.2f} "
                     f"between methods", fontsize=10)
        ax.grid(alpha=0.22, lw=0.6)
        ax.set_axisbelow(True)
    handles = [Line2D([], [], marker="o", ls="", ms=8, mfc=C_BINDER, mec="k",
                      label="binder"),
               Line2D([], [], marker="o", ls="", ms=6, mfc=C_DECOY,
                      mec="white", label="non-binder")]
    fig.suptitle("Do the two methods see the same thing? "
                 "(dashed lines = per-panel medians; boxes = binders / designs "
                 "per quadrant)", y=1.06, fontsize=11.5)
    fig.legend(handles=handles, loc="upper center", ncol=2, fontsize=9,
               frameon=False, bbox_to_anchor=(0.5, 1.005))
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(FIG, "fig4_method_agreement.png")
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("wrote", os.path.relpath(out, pp.ROOT))


# ------------------------------------------------------------------ figure 5
def fig5(store):
    """Screening yield: order the top X% of the library by each score and see
    what fraction of the true binders you captured (recall) and what fraction
    of your picks were real (precision)."""
    strategies = [(f"{E} ipTM", "#4C72B0", "-"),
                  (f"{P} ipTM", "#DD8452", "-"),
                  ("mean-z(ipTM)", "#55A868", "-"),
                  ("min-z(ipTM)", "#55A868", "--"),
                  (f"{E} interface pLDDT", "#4C72B0", ":")]
    fig, axes = plt.subplots(2, 3, figsize=(14, 7.2), sharex=True)
    fracs = np.linspace(0.02, 1.0, 60)
    for j, task in enumerate(TASK_ORDER):
        df, allp = get_scores(store, task)
        y = df["label"].values.astype(int)
        n, prev = len(y), y.mean()
        for key, c, ls in strategies:
            if key not in allp:
                continue
            order = np.argsort(-allp[key], kind="stable")
            ys = y[order]
            ks = np.maximum(1, np.round(fracs * n).astype(int))
            rec = np.array([ys[:k].sum() / y.sum() for k in ks])
            prc = np.array([ys[:k].mean() for k in ks])
            axes[0, j].plot(fracs, rec, color=c, ls=ls, lw=1.8, label=key)
            axes[1, j].plot(fracs, prc / prev, color=c, ls=ls, lw=1.8, label=key)
        axes[0, j].plot(fracs, fracs, color="k", ls="--", lw=0.9, alpha=0.5,
                        label="random")
        axes[1, j].axhline(1.0, color="k", ls="--", lw=0.9, alpha=0.5,
                           label="random")
        axes[0, j].set_title(f"{TASK_LABEL[task]}\nprevalence {prev:.0%} "
                             f"({int(y.sum())}/{n})", fontsize=10)
        axes[1, j].set_xlabel("fraction of the 95-design library ordered")
        for r in (0, 1):
            axes[r, j].grid(alpha=0.25, lw=0.6)
            axes[r, j].set_axisbelow(True)
    axes[0, 0].set_ylabel("recall\n(fraction of true binders captured)")
    axes[1, 0].set_ylabel("enrichment factor\n(precision / prevalence)")
    axes[0, 0].legend(fontsize=7.5, loc="lower right")
    fig.suptitle("Screening yield: what you get if you only order the top-ranked "
                 "designs", y=0.99, fontsize=12)
    fig.tight_layout()
    out = os.path.join(FIG, "fig5_screening_yield.png")
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("wrote", os.path.relpath(out, pp.ROOT))


# ------------------------------------------------------------------ figure 6
def fig6(budget_df):
    """Budget-matched strategies: for a fixed number of designs you can afford
    to express, which way of choosing them finds the most real binders?"""
    strat_order = [f"{E} top-k", f"{P} top-k", "fused mean-z top-k",
                   "union: top-k/2 from each", "AND gate: top-k in both"]
    colors = ["#4C72B0", "#DD8452", "#55A868", "#8172B2", "#937860"]
    budgets = sorted(budget_df.budget_k.unique())
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharey=False)
    w = 0.8 / len(strat_order)
    for ax, task in zip(axes, TASK_ORDER):
        sub = budget_df[budget_df.task == task]
        npos = int(sub["n_pos_total"].iloc[0])
        x = np.arange(len(budgets))
        for si, (strat, c) in enumerate(zip(strat_order, colors)):
            vals, picked = [], []
            for k in budgets:
                r = sub[(sub.budget_k == k) & (sub.strategy == strat)]
                vals.append(float(r["n_binders_found"].iloc[0]) if len(r) else 0.0)
                picked.append(int(r["n_picked"].iloc[0]) if len(r) else 0)
            pos = x - 0.4 + w * (si + 0.5)
            ax.bar(pos, vals, width=w * 0.92, color=c, alpha=0.9,
                   label=strat, edgecolor="white", linewidth=0.5)
            for xi, (v, nk) in enumerate(zip(vals, picked)):
                if strat.startswith("AND") and nk != budgets[xi]:
                    ax.text(pos[xi], v + 0.12, f"n={nk}", ha="center",
                            fontsize=6.2, color="#555555", rotation=90)
        ax.axhline(npos, color="k", ls="--", lw=0.9, alpha=0.6)
        ax.text(len(budgets) - 0.55, npos, f" all {npos}", va="bottom",
                ha="right", fontsize=7.5, color="k")
        ax.set_xticks(x)
        ax.set_xticklabels([str(k) for k in budgets])
        ax.set_xlabel("budget k  (designs expressed and tested)")
        ax.set_ylabel("true binders found")
        ax.set_title(f"{TASK_LABEL[task]}\n{npos} binders in "
                     f"{int(sub['n_designs'].iloc[0])} designs", fontsize=10)
        ax.grid(axis="y", alpha=0.25, lw=0.6)
        ax.set_axisbelow(True)
    axes[0].legend(fontsize=7.5, loc="upper left")
    fig.suptitle("Spending a fixed screening budget: one method, both fused, "
                 "or the budget split between them", y=1.02, fontsize=12)
    fig.tight_layout()
    out = os.path.join(FIG, "fig6_budget_matched.png")
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("wrote", os.path.relpath(out, pp.ROOT))


def main():
    wide, tasks, store, auc_df, h2h_df, gate_df, comp_df, budget_df = pp.main()
    print()
    fig1(store)
    fig2(store)
    fig3(auc_df)
    fig4(wide, tasks)
    fig5(store)
    fig6(budget_df)


if __name__ == "__main__":
    main()
