#!/usr/bin/env python3
"""
Predictive power of ESMFold2 and Protenix-v2 reprediction -- alone and together.

Ground truth: BLI binding at 1 uM (PD1_duo_DesignSummary.csv -> Classification).
Three retrospective screening tasks:

  human           : does the design bind hPD-1?   (binder vs decoy)
  mouse           : does the design bind mPD-1?   (binder vs decoy)
  cross_reactive  : does the design bind BOTH?    (cross-reactive vs single/none)
                    one score per design, taken as the worst of its two species
                    scores -- a design is only cross-reactive if its weaker side
                    still holds up.

For every candidate predictor we report
  AUC                 with a stratified bootstrap 95% CI
  AUC_0.1             partial AUC restricted to FPR <= 0.1, normalised
  P@10% / P@20%, EF   precision and enrichment in the top-ranked fraction
                      (what actually matters when you can only order 10 designs)
and, for the head-to-head questions, a DeLong test on the paired DeltaAUC.

"Together" covers three genuinely different ways to combine the two methods:
  score fusion        mean-z / min / rank-mean of the two ipTMs, and a
                      cross-validated logistic model
  structural consensus  do the two models place the binder on the same epitope
                      in the same pose (cross-method RMSD / epitope Jaccard)?
                      This feature does not exist for either method alone.
  gated selection     require both methods to pass -- measured as precision at
                      a matched number of picks.

The natural-ligand control (HuPDL1-Fc) is excluded from all tasks and reported
separately: it is a 227 aa PDB-deposited complex, not a de novo design.
"""

import os
import json
import itertools
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_curve, roc_auc_score, auc as sk_auc

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "Analysis", "data")
FIG = os.path.join(ROOT, "Analysis", "figures")
os.makedirs(FIG, exist_ok=True)

RNG = np.random.default_rng(0)
N_BOOT = 4000
PARTIAL_FPR_MAX = 0.10
TOP_FRACS = [0.10, 0.20]
CV_REPEATS = 20
CV_FOLDS = 5

# build_dataset.py already removes the natural-ligand positive control from
# every table it writes. These guards are belt-and-braces: they keep the
# analysis correct if it is ever pointed at a table that still contains it.
POSITIVE_CONTROL = "HuPDL1-Fc"

# Classification -> label per species. None = excluded (ambiguous / no data)
CLASS_TO_LABEL = {
    "cross-reactive":                {"human": 1,    "mouse": 1},
    "human-specific":                {"human": 1,    "mouse": 0},
    "mouse-specific":                {"human": 0,    "mouse": 1},
    "non-binding":                   {"human": 0,    "mouse": 0},
    "non-binding (no mupd-1 data)":  {"human": 0,    "mouse": None},
    "likely nonspecific/sticky":     {"human": 0,    "mouse": 0},
    "borderline (weak human trend)": {"human": None, "mouse": 0},
    "borderline (weak mouse trend)": {"human": 0,    "mouse": None},
}

# ---------------------------------------------------------------------------
# Predictor registry.  kind:
#   'col'   single column, sign = +1 if higher is better
#   'fuse'  parameter-free combination of columns (no fitting, no labels)
#   'cv'    cross-validated logistic regression over several columns
# ---------------------------------------------------------------------------
E, P = "ESMFold2", "Protenix-v2"

SINGLE_METRICS = [
    ("ipTM",               +1, "ipTM"),
    ("pTM",                +1, "pTM"),
    ("pLDDT",              +1, "pLDDT"),
    ("pAE",                -1, "pAE"),
    ("i_pAE",              -1, "i-pAE"),
    ("iface_pLDDT_binder", +1, "interface pLDDT"),
    ("n_contacts",         +1, "interface contacts"),
    ("native_epitope_jaccard", +1, "PD-L1-epitope overlap"),
    ("xs_epitope_jaccard", +1, "h/m epitope conservation"),
]

CV_FEATURES = {
    f"{E} composite (CV)":  [f"{E}_{c}" for c in
                             ["ipTM", "pTM", "pLDDT", "pAE", "iface_pLDDT_binder"]],
    f"{P} composite (CV)":  [f"{P}_{c}" for c in
                             ["ipTM", "pTM", "pLDDT", "pAE", "i_pAE", "iface_pLDDT_binder"]],
    "Both ipTM (CV)":       [f"{E}_ipTM", f"{P}_ipTM"],
    "Both ipTM + consensus (CV)": [f"{E}_ipTM", f"{P}_ipTM",
                                   "xm_binder_rmsd", "xm_epitope_jaccard"],
    "All metrics, both methods (CV)":
        [f"{E}_{c}" for c in ["ipTM", "pTM", "pLDDT", "pAE", "iface_pLDDT_binder"]] +
        [f"{P}_{c}" for c in ["ipTM", "pTM", "pLDDT", "pAE", "i_pAE", "iface_pLDDT_binder"]],
}


# ============================================================== ROC utilities
def partial_auc(y, s):
    fpr, tpr, _ = roc_curve(y, s)
    if PARTIAL_FPR_MAX >= fpr.max():
        return sk_auc(fpr, tpr) * min(PARTIAL_FPR_MAX, fpr.max()) / PARTIAL_FPR_MAX
    within = fpr <= PARTIAL_FPR_MAX
    fx, ty = fpr[within].tolist(), tpr[within].tolist()
    fx.append(PARTIAL_FPR_MAX)
    ty.append(float(np.interp(PARTIAL_FPR_MAX, fpr, tpr)))
    trapz = getattr(np, "trapezoid", None) or np.trapz
    return trapz(ty, fx) / PARTIAL_FPR_MAX


def topk_metrics(y, s):
    n = len(y)
    prev = y.mean()
    order = np.argsort(-s, kind="stable")
    ys = y[order]
    out = {"prevalence": prev, "n": n, "n_pos": int(y.sum())}
    for f in TOP_FRACS:
        k = max(1, int(round(f * n)))
        pk = ys[:k].mean()
        out[f"P@{int(f*100)}%"] = pk
        out[f"EF@{int(f*100)}%"] = pk / prev if prev > 0 else np.nan
        out[f"recall@{int(f*100)}%"] = ys[:k].sum() / y.sum() if y.sum() else np.nan
        out[f"k@{int(f*100)}%"] = k
    return out


def bootstrap_aucs(y, s, n_boot=N_BOOT, rng=None):
    """Stratified bootstrap AUCs, fully vectorised. Positives and negatives are
    resampled separately so every replicate keeps both classes and the class
    balance fixed. AUC comes from the Mann-Whitney rank identity, with midranks
    handling the ties that resampling inevitably creates."""
    rng = rng or RNG
    pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
    m, n = len(pos), len(neg)
    if m < 2 or n < 2:
        return np.array([])
    idx = np.concatenate([rng.choice(pos, (n_boot, m), replace=True),
                          rng.choice(neg, (n_boot, n), replace=True)], axis=1)
    R = stats.rankdata(s[idx], axis=1)
    return (R[:, :m].sum(axis=1) - m * (m + 1) / 2.0) / (m * n)


def bootstrap_auc_ci(y, s, n_boot=N_BOOT, rng=None):
    a = bootstrap_aucs(y, s, n_boot, rng)
    if a.size == 0:
        return np.nan, np.nan
    return float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))


# ------------------------------------------------------------ DeLong (fast)
def _midrank(x):
    J = np.argsort(x)
    Z = x[J]
    N = len(x)
    T = np.zeros(N)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1)
        i = j
    T2 = np.empty(N)
    T2[J] = T + 1
    return T2


def delong_cov(scores, y):
    """Fast DeLong (Sun & Xu 2014). scores: k x n. Returns (aucs, k x k cov)."""
    order = np.argsort(-y, kind="stable")          # positives first
    s = np.asarray(scores, float)[:, order]
    m = int(y.sum())
    n = s.shape[1] - m
    k = s.shape[0]
    tx = np.vstack([_midrank(s[r, :m]) for r in range(k)])
    ty = np.vstack([_midrank(s[r, m:]) for r in range(k)])
    tz = np.vstack([_midrank(s[r, :]) for r in range(k)])
    aucs = tz[:, :m].sum(axis=1) / m / n - (m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx) / n
    v10 = 1.0 - (tz[:, m:] - ty) / m
    sx = np.atleast_2d(np.cov(v01))
    sy = np.atleast_2d(np.cov(v10))
    return aucs, sx / m + sy / n


def delong_test(y, s1, s2):
    """Two-sided p-value for AUC(s1) == AUC(s2) on the same samples."""
    aucs, cov = delong_cov(np.vstack([s1, s2]), y)
    l = np.array([1.0, -1.0])
    var = float(l @ cov @ l)
    if var <= 0:
        return aucs[0] - aucs[1], np.nan, np.nan
    z = (aucs[0] - aucs[1]) / np.sqrt(var)
    return float(aucs[0] - aucs[1]), float(z), float(2 * stats.norm.sf(abs(z)))


# ================================================================ data layer
N_DESIGNS = None   # set by load_tasks() from the merged tables


def load_tasks():
    global N_DESIGNS
    wide = {sp: pd.read_csv(os.path.join(DATA, f"merged_{sp}.csv"))
            for sp in ("human", "mouse")}
    N_DESIGNS = wide["human"]["name"].nunique()
    for sp, d in wide.items():
        d["Classification"] = d["Classification"].astype(str).str.strip()
        d["label"] = [CLASS_TO_LABEL[c.lower()][sp] for c in d["Classification"]]
        d.loc[d["name"] == POSITIVE_CONTROL, "label"] = np.nan

    # ---- per-species tasks
    tasks = {}
    for sp in ("human", "mouse"):
        d = wide[sp]
        keep = d["label"].notna() & (d["name"] != POSITIVE_CONTROL)
        tasks[sp] = d.loc[keep].reset_index(drop=True).copy()
        tasks[sp]["label"] = tasks[sp]["label"].astype(int)

    # ---- cross-reactivity task: one row per design, worst-of-two-species
    h, m = wide["human"], wide["mouse"]
    common = [c for c in h.columns if c in m.columns and
              (h[c].dtype.kind in "fi" and m[c].dtype.kind in "fi")]
    cr = h[["name", "VisterraName", "Classification"]].copy()
    signs = {}
    for c in common:
        if c in ("label",):
            continue
        sign = -1 if ("pAE" in c or "rmsd" in c) else +1
        signs[c] = sign
        # worst-case = min for higher-is-better, max for lower-is-better
        cr[c] = (np.minimum(h[c].values, m[c].values) if sign > 0
                 else np.maximum(h[c].values, m[c].values))
    cr["label"] = (cr["Classification"].str.lower() == "cross-reactive").astype(int)
    cr = cr[cr["name"] != POSITIVE_CONTROL].reset_index(drop=True)
    tasks["cross_reactive"] = cr
    return wide, tasks


# ============================================================ score builders
def cv_scores(df, cols, label_col="label", n_repeats=CV_REPEATS):
    """Out-of-fold logistic-regression scores. Returns (mean OOF score,
    list of per-repeat OOF AUCs). Standardisation is fitted inside each fold,
    so no information leaks from held-out designs."""
    X = df[cols].values.astype(float)
    y = df[label_col].values.astype(int)
    ok = np.isfinite(X).all(axis=1)
    if not ok.all():
        return None, None, ok
    oof = np.zeros((n_repeats, len(y)))
    cv = RepeatedStratifiedKFold(n_splits=CV_FOLDS, n_repeats=n_repeats,
                                 random_state=0)
    aucs = []
    for i, (tr, te) in enumerate(cv.split(X, y)):
        rep = i // CV_FOLDS
        mdl = make_pipeline(StandardScaler(),
                            LogisticRegression(C=1.0, max_iter=5000))
        mdl.fit(X[tr], y[tr])
        oof[rep, te] = mdl.predict_proba(X[te])[:, 1]
    for rep in range(n_repeats):
        aucs.append(roc_auc_score(y, oof[rep]))
    return oof.mean(axis=0), np.array(aucs), ok


def zscore(v):
    v = np.asarray(v, float)
    sd = np.nanstd(v)
    return (v - np.nanmean(v)) / (sd if sd > 0 else 1.0)


def build_predictors(df):
    """Return {label: (score_vector, group, note)} for one task table."""
    out = {}

    def add(lbl, vec, group, note=""):
        v = np.asarray(vec, float)
        if np.all(np.isnan(v)) or np.nanstd(v) == 0:
            return
        out[lbl] = (v, group, note)

    # --- singles, per method
    for meth, group in ((E, f"{E} alone"), (P, f"{P} alone")):
        for col, sign, pretty in SINGLE_METRICS:
            c = f"{meth}_{col}"
            if c in df.columns and df[c].notna().any():
                add(f"{meth} {pretty}", sign * df[c].values, group)

    # --- parameter-free fusion of the two ipTMs (uses no labels)
    e, p = f"{E}_ipTM", f"{P}_ipTM"
    if e in df.columns and p in df.columns:
        add("mean-z(ipTM)", 0.5 * (zscore(df[e]) + zscore(df[p])),
            "Together (fusion)", "unweighted z-average")
        add("min-z(ipTM)", np.minimum(zscore(df[e]), zscore(df[p])),
            "Together (fusion)", "both must be high (AND-like)")
        add("max-z(ipTM)", np.maximum(zscore(df[e]), zscore(df[p])),
            "Together (fusion)", "either may be high (OR-like)")
        r = 0.5 * (stats.rankdata(df[e]) + stats.rankdata(df[p]))
        add("rank-mean(ipTM)", r, "Together (fusion)", "scale-free")

    # --- cross-method structural consensus (exists only jointly)
    if "xm_binder_rmsd" in df.columns:
        add("cross-method pose RMSD", -df["xm_binder_rmsd"].values,
            "Together (consensus)", "lower = the two models agree")
    if "xm_epitope_jaccard" in df.columns:
        add("cross-method epitope overlap", df["xm_epitope_jaccard"].values,
            "Together (consensus)", "")
    if "xm_paratope_jaccard" in df.columns:
        add("cross-method paratope overlap", df["xm_paratope_jaccard"].values,
            "Together (consensus)", "")
    if "xm_binder_fold_rmsd" in df.columns:
        add("cross-method fold RMSD", -df["xm_binder_fold_rmsd"].values,
            "Together (consensus)", "same fold, ignoring pose")

    # --- baselines
    if "length" in df.columns:
        add("binder length (baseline)", df["length"].values, "Baseline",
            "confound check")

    return out


def build_cv_predictors(df):
    """Cross-validated multivariate models. Returned separately because their
    AUC is an out-of-fold estimate, not an apparent one."""
    out = {}
    for lbl, cols in CV_FEATURES.items():
        cols_ok = [c for c in cols if c in df.columns and df[c].notna().all()]
        if len(cols_ok) < 2:
            continue
        score, aucs, ok = cv_scores(df, cols_ok)
        if score is None:
            continue
        group = ("Together (learned)" if any(c.startswith(E) for c in cols_ok)
                 and any(c.startswith(P) for c in cols_ok)
                 else (f"{E} alone" if cols_ok[0].startswith(E) else f"{P} alone"))
        out[lbl] = (score, group, aucs, cols_ok)
    return out


# ==================================================================== report
def evaluate(task_name, df):
    y = df["label"].values.astype(int)
    rows = []

    preds = build_predictors(df)
    for lbl, (s, group, note) in preds.items():
        ok = np.isfinite(s)
        yy, ss = y[ok], s[ok]
        if len(np.unique(yy)) < 2:
            continue
        a = roc_auc_score(yy, ss)
        lo, hi = bootstrap_auc_ci(yy, ss)
        rows.append({"task": task_name, "predictor": lbl, "group": group,
                     "AUC": a, "AUC_lo": lo, "AUC_hi": hi,
                     "AUC_0.1": partial_auc(yy, ss), "estimate": "apparent",
                     "note": note, **topk_metrics(yy, ss)})

    cvp = build_cv_predictors(df)
    for lbl, (s, group, aucs, cols_ok) in cvp.items():
        a = roc_auc_score(y, s)
        rows.append({"task": task_name, "predictor": lbl, "group": group,
                     "AUC": float(np.mean(aucs)),
                     "AUC_lo": float(np.percentile(aucs, 2.5)),
                     "AUC_hi": float(np.percentile(aucs, 97.5)),
                     "AUC_0.1": partial_auc(y, s), "estimate": "out-of-fold",
                     "note": f"{len(cols_ok)} features, {CV_REPEATS}x{CV_FOLDS}-fold CV "
                             f"(CI = spread over repeats)",
                     **topk_metrics(y, s)})

    return pd.DataFrame(rows), preds, cvp


def head_to_head(task_name, y, preds, cvp):
    """DeLong tests for the questions that actually matter."""
    allp = {k: v[0] for k, v in preds.items()}
    allp.update({k: v[0] for k, v in cvp.items()})
    e_best, p_best = f"{E} ipTM", f"{P} ipTM"
    pairs = [
        (p_best, e_best,                     "Protenix-v2 ipTM vs ESMFold2 ipTM"),
        ("mean-z(ipTM)", e_best,             "fused ipTM vs ESMFold2 alone"),
        ("mean-z(ipTM)", p_best,             "fused ipTM vs Protenix-v2 alone"),
        ("Both ipTM (CV)", e_best,           "learned fusion vs ESMFold2 alone"),
        ("Both ipTM + consensus (CV)", "Both ipTM (CV)",
                                             "adding structural consensus to fusion"),
        ("All metrics, both methods (CV)", f"{E} composite (CV)",
                                             "both methods vs ESMFold2 metrics only"),
        (f"{E} interface pLDDT", e_best,     "interface pLDDT vs ipTM (ESMFold2)"),
        (f"{P} interface pLDDT", p_best,     "interface pLDDT vs ipTM (Protenix-v2)"),
    ]
    rows = []
    for a, b, q in pairs:
        if a not in allp or b not in allp:
            continue
        sa, sb = allp[a], allp[b]
        ok = np.isfinite(sa) & np.isfinite(sb)
        d, z, pv = delong_test(y[ok], sa[ok], sb[ok])
        rows.append({"task": task_name, "question": q, "A": a, "B": b,
                     "AUC_A": roc_auc_score(y[ok], sa[ok]),
                     "AUC_B": roc_auc_score(y[ok], sb[ok]),
                     "dAUC": d, "z": z, "p_DeLong": pv})
    return pd.DataFrame(rows)


def gated_selection(task_name, df):
    """'Together' as a filter, not a score: keep designs that clear the same
    percentile in BOTH methods, then compare the precision/recall of that pick
    against each method's own top-N and against the fused score at a matched N
    -- i.e. at an equal ordering budget. Tidy long format."""
    y = df["label"].values.astype(int)
    e, p = df[f"{E}_ipTM"].values, df[f"{P}_ipTM"].values
    fused = 0.5 * (zscore(e) + zscore(p))
    rows = []
    for q in [0.5, 0.6, 0.7, 0.75, 0.8, 0.9]:
        sel = (e >= np.quantile(e, q)) & (p >= np.quantile(p, q))
        k = int(sel.sum())
        if k == 0:
            continue
        strategies = {
            "both methods (AND gate)": sel,
            f"{E} top-N": np.isin(np.arange(len(y)),
                                  np.argsort(-e, kind="stable")[:k]),
            f"{P} top-N": np.isin(np.arange(len(y)),
                                  np.argsort(-p, kind="stable")[:k]),
            "fused mean-z top-N": np.isin(np.arange(len(y)),
                                          np.argsort(-fused, kind="stable")[:k]),
        }
        for name, mask in strategies.items():
            rows.append({"task": task_name, "percentile": q, "n_picked": k,
                         "frac_library": k / len(y), "strategy": name,
                         "precision": y[mask].mean(),
                         "recall": y[mask].sum() / y.sum(),
                         "n_true_binders": int(y[mask].sum()),
                         "prevalence": y.mean(),
                         "EF": y[mask].mean() / y.mean() if y.mean() else np.nan})
    return pd.DataFrame(rows)


def budget_matched_strategies(task_name, df, budgets=(10, 15, 20, 25, 30)):
    """The fair way to ask whether two methods beat one: fix the number of
    designs you can afford to express and test (the budget k), then compare

      ESMFold2 top-k            one method, whole budget
      Protenix-v2 top-k         the other method, whole budget
      fused mean-z top-k        one ranked list built from both scores
      union, k/2 each           split the budget: top-k/2 from each method,
                                take the union (deduplicated, then topped up
                                from the fused ranking so every strategy
                                really does spend the same k)
      AND gate, top-k of both   only designs each method ranks in its own top-k

    Recall (how many of the true binders you captured) is the number that
    matters for a discovery campaign; precision is what it costs you.
    """
    y = df["label"].values.astype(int)
    e, p = df[f"{E}_ipTM"].values, df[f"{P}_ipTM"].values
    fused = 0.5 * (zscore(e) + zscore(p))
    oe = np.argsort(-e, kind="stable")
    op = np.argsort(-p, kind="stable")
    of = np.argsort(-fused, kind="stable")
    n, npos = len(y), int(y.sum())
    rows = []
    for k in budgets:
        if k > n:
            continue
        half = int(np.ceil(k / 2))
        union = list(dict.fromkeys(list(oe[:half]) + list(op[:half])))
        for extra in of:                      # top the union back up to k
            if len(union) >= k:
                break
            if extra not in union:
                union.append(extra)
        both = [i for i in oe[:k] if i in set(op[:k])]
        picks = {
            f"{E} top-k": list(oe[:k]),
            f"{P} top-k": list(op[:k]),
            "fused mean-z top-k": list(of[:k]),
            "union: top-k/2 from each": union[:k],
            "AND gate: top-k in both": both,
        }
        for name, idx in picks.items():
            idx = np.array(idx, dtype=int)
            if len(idx) == 0:
                continue
            rows.append({"task": task_name, "budget_k": k, "strategy": name,
                         "n_picked": len(idx),
                         "n_binders_found": int(y[idx].sum()),
                         "recall": y[idx].sum() / npos,
                         "precision": y[idx].mean(),
                         "EF": y[idx].mean() / y.mean(),
                         "n_pos_total": npos, "n_designs": n})
    return pd.DataFrame(rows)


def power_note(h2h_df):
    """How big a difference could this dataset even detect? The DeLong standard
    error is recoverable from each test (SE = dAUC / z), so the minimum
    detectable DeltaAUC at 80% power / alpha=0.05 two-sided is 2.80 * SE.
    Anything smaller than that is not evidence of equivalence -- just silence."""
    rows = []
    for _, r in h2h_df.iterrows():
        if not np.isfinite(r["z"]) or r["z"] == 0:
            continue
        se = abs(r["dAUC"] / r["z"])
        rows.append({"task": r["task"], "question": r["question"],
                     "dAUC": r["dAUC"], "SE_DeLong": se,
                     "min_detectable_dAUC_80pct": 2.80 * se})
    return pd.DataFrame(rows)


def complementarity(tasks, wide):
    rows = []
    for sp in ("human", "mouse"):
        d = wide[sp]
        for col in ["ipTM", "pTM", "pLDDT", "pAE", "iface_pLDDT_binder"]:
            a, b = d[f"{E}_{col}"], d[f"{P}_{col}"]
            ok = a.notna() & b.notna() & (d["name"] != POSITIVE_CONTROL)
            rows.append({"species": sp, "metric": col,
                         "pearson_r": stats.pearsonr(a[ok], b[ok])[0],
                         "spearman_rho": stats.spearmanr(a[ok], b[ok])[0],
                         "n": int(ok.sum())})
        ok = (d["name"] != POSITIVE_CONTROL)
        rows.append({"species": sp, "metric": "xm_binder_rmsd (median A)",
                     "pearson_r": np.nan,
                     "spearman_rho": float(d.loc[ok, "xm_binder_rmsd"].median()),
                     "n": int(ok.sum())})
    return pd.DataFrame(rows)


def main():
    wide, tasks = load_tasks()

    print("=" * 78)
    print("TASK DEFINITIONS")
    print("=" * 78)
    for t, d in tasks.items():
        y = d["label"].values
        print(f"  {t:15s} n={len(d):3d}  positives={int(y.sum()):3d} "
              f"({y.mean():.1%})   excluded={N_DESIGNS - len(d)}/"
              f"{N_DESIGNS} as ambiguous or no-data for this task")
    print("  (the HuPDL1-Fc natural-ligand positive control is not in the "
          "dataset at all -- removed by build_dataset.py)")

    all_auc, all_h2h, all_gate, all_budget = [], [], [], []
    store = {}
    for t, d in tasks.items():
        res, preds, cvp = evaluate(t, d)
        all_auc.append(res)
        all_h2h.append(head_to_head(t, d["label"].values.astype(int), preds, cvp))
        all_gate.append(gated_selection(t, d))
        all_budget.append(budget_matched_strategies(t, d))
        store[t] = (d, preds, cvp)
    budget_df = pd.concat(all_budget, ignore_index=True)

    auc_df = pd.concat(all_auc, ignore_index=True)
    h2h_df = pd.concat(all_h2h, ignore_index=True)
    gate_df = pd.concat(all_gate, ignore_index=True)
    comp_df = complementarity(tasks, wide)

    auc_df.to_csv(os.path.join(DATA, "auc_table.csv"), index=False)
    h2h_df.to_csv(os.path.join(DATA, "head_to_head.csv"), index=False)
    gate_df.to_csv(os.path.join(DATA, "gated_selection.csv"), index=False)
    comp_df.to_csv(os.path.join(DATA, "method_correlation.csv"), index=False)
    budget_df.to_csv(os.path.join(DATA, "budget_matched.csv"), index=False)
    rescue_table(tasks, os.path.join(DATA, "binder_ranks.csv"))

    show = ["predictor", "group", "estimate", "AUC", "AUC_lo", "AUC_hi",
            "AUC_0.1", "P@10%", "EF@10%", "recall@10%", "n_pos", "n"]
    for t in tasks:
        sub = auc_df[auc_df.task == t].sort_values("AUC", ascending=False)
        print("\n" + "=" * 78)
        print(f"AUC RANKING -- task: {t}   (prevalence "
              f"{sub['prevalence'].iloc[0]:.1%})")
        print("=" * 78)
        print(sub[show].to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    print("\n" + "=" * 78)
    print("HEAD-TO-HEAD (DeLong test on paired DeltaAUC)")
    print("=" * 78)
    print(h2h_df.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    pw = power_note(h2h_df)
    pw.to_csv(os.path.join(DATA, "power.csv"), index=False)
    print("\n" + "=" * 78)
    print("WHAT COULD THIS DATASET EVEN DETECT? (min detectable DeltaAUC, "
          "80% power, alpha=0.05)")
    print("=" * 78)
    print(pw.groupby("task")["min_detectable_dAUC_80pct"]
          .agg(["min", "median", "max"])
          .to_string(float_format=lambda x: f"{x:.3f}"))

    print("\n" + "=" * 78)
    print("GATED SELECTION -- require BOTH methods above the same percentile,")
    print("compared with each method's own top-N at a matched N")
    print("=" * 78)
    print(gate_df.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    print("\n" + "=" * 78)
    print("BUDGET-MATCHED STRATEGIES -- same number of designs expressed, "
          "different way of choosing them")
    print("=" * 78)
    for t in tasks:
        sub = budget_df[budget_df.task == t]
        piv = sub.pivot_table(index="strategy", columns="budget_k",
                              values="n_binders_found")
        print(f"\n  {t}  (true binders found, out of "
              f"{int(sub['n_pos_total'].iloc[0])}; columns = budget k)")
        print(piv.to_string(float_format=lambda x: f"{x:.0f}"))

    print("\n" + "=" * 78)
    print("HOW CORRELATED ARE THE TWO METHODS? (designs only, control excluded)")
    print("=" * 78)
    print(comp_df.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    ctrl_path = os.path.join(DATA, "positive_control.csv")
    if os.path.isfile(ctrl_path):
        ctrl = pd.read_csv(ctrl_path)
        print("\n" + "=" * 78)
        print("POSITIVE CONTROL (HuPDL1-Fc, natural ligand) -- REMOVED from the "
              "dataset, shown for reference only")
        print("=" * 78)
        print(ctrl[["species", "method", "i_pTM", "pLDDT", "iface_pLDDT_binder",
                    "xm_binder_rmsd"]].to_string(
              index=False, float_format=lambda x: f"{x:.3f}"))

    with open(os.path.join(DATA, "summary.json"), "w") as fh:
        json.dump({"tasks": {t: {"n": int(len(d)),
                                 "n_pos": int(d['label'].sum())}
                             for t, d in tasks.items()}}, fh, indent=2)
    print(f"\nWrote tables to Analysis/data/: auc_table.csv, head_to_head.csv, "
          f"gated_selection.csv, method_correlation.csv")
    return wide, tasks, store, auc_df, h2h_df, gate_df, comp_df, budget_df


if __name__ == "__main__":
    main()


# ---------------------------------------------------------------------------
# Design-level complementarity: are there true binders that ONE method ranks
# highly and the other misses? That is the only way fusion can ever help.
# ---------------------------------------------------------------------------
def rescue_table(tasks, out_csv=None):
    rows = []
    for task, df in tasks.items():
        y = df["label"].values.astype(int)
        e, p = df[f"{E}_ipTM"].values, df[f"{P}_ipTM"].values
        n = len(y)
        # percentile rank, 1.0 = best
        re_ = stats.rankdata(e) / n
        rp = stats.rankdata(p) / n
        for i in np.where(y == 1)[0]:
            rows.append({"task": task, "name": df["name"].iloc[i],
                         "VisterraName": df["VisterraName"].iloc[i],
                         "Classification": df["Classification"].iloc[i],
                         f"{E}_ipTM": e[i], f"{E}_pctile": re_[i],
                         f"{P}_ipTM": p[i], f"{P}_pctile": rp[i],
                         "pctile_gap": re_[i] - rp[i]})
    t = pd.DataFrame(rows)
    if out_csv:
        t.to_csv(out_csv, index=False)
    return t
