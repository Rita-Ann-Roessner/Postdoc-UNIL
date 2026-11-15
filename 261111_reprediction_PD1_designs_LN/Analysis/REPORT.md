# Predictive power of ESMFold2 and Protenix-v2 reprediction — alone and together

PD1_duo campaign, 95 de novo designs, BLI at 1 µM as ground truth. The HuPDL1-Fc
natural-ligand control is removed from the dataset entirely.

## What was benchmarked

| | |
|---|---|
| Designs | 95 (the HuPDL1-Fc natural-ligand control is removed from every table) |
| Predictions | ESMFold2 and Protenix-v2, each against hPD-1 and mPD-1 → 380 complexes |
| Ground truth | `PD1_duo_DesignSummary.csv` → `Classification` from BLI at 1 µM |
| Tasks | hPD-1 binding (16 binders / 94), mPD-1 binding (12 / 93), cross-reactivity (9 / 95) |
| Scores tested | 9 pipeline/derived metrics per method, 4 parameter-free fusions, 5 cross-validated models, 4 cross-method consensus features, 1 confound baseline |

Everything is recomputed from the raw outputs — including per-residue interface pLDDT,
interface contacts, epitope identity and cross-method pose RMSD, parsed straight out of the
384 predicted complex PDBs.

## Headline

| Score | hPD-1 | mPD-1 | cross-reactivity |
|---|---|---|---|
| **ESMFold2 ipTM** | 0.758 [0.64–0.86] | 0.749 [0.63–0.85] | **0.819 [0.68–0.93]** |
| **ESMFold2 interface pLDDT** | **0.787 [0.67–0.89]** | 0.706 [0.54–0.85] | 0.758 [0.59–0.90] |
| **Protenix-v2 ipTM** | 0.640 [0.49–0.77] | 0.746 [0.61–0.87] | 0.661 [0.46–0.83] |
| mean-z fusion of both ipTMs | 0.730 [0.59–0.85] | **0.776 [0.67–0.87]** | 0.771 [0.58–0.92] |
| logistic model on both ipTMs (out-of-fold) | 0.726 | 0.745 | 0.762 |
| cross-method pose agreement | 0.593 | 0.607 | 0.534 |
| binder length (confound check) | 0.472 | 0.496 | 0.444 |

AUC with stratified-bootstrap 95% CI. Full table: `Analysis/data/auc_table.csv`.

## Four findings

**1. Both methods work, weakly.** Every ipTM-type score separates binders from non-binders
better than chance (p = 0.0006–0.06, one-sided Mann-Whitney). But at a realistic screening
budget the enrichment is only ~2–3× over random: ESMFold2's top 10 hPD-1 picks contain 3 real
binders out of 16, its top 25 contain 9. This is a useful triage filter, not a decision procedure.

**2. ESMFold2 ≥ Protenix-v2, and the gap is on human PD-1.** ESMFold2 wins hPD-1 by
ΔAUC 0.118 (DeLong p = 0.060) and cross-reactivity by ΔAUC 0.158 (p = 0.022). On mPD-1 they
are indistinguishable (ΔAUC 0.003, p = 0.97). Protenix-v2's ipTM is also compressed into a
narrow high band (0.84–0.93 for most designs), which is what costs it discrimination.

**3. Recomputed interface pLDDT is at least as good as ipTM — and free.** Averaging
per-residue pLDDT over just the interface residues gives the best single hPD-1 predictor in
the whole benchmark (0.787). It is not in the pipeline's metrics table; it comes from the
B-factor column of PDBs you already have.

**4. Combining the two methods does not help, by any of three routes.**
- *Score fusion*: mean-z, min-z, max-z and rank-mean all land inside the CI of the better
  single method. Best case (mouse max-z, 0.785 vs 0.749) is ΔAUC 0.036 against a minimum
  detectable difference of ~0.21 for that task.
- *Learned fusion*: out-of-fold logistic regression on both ipTMs gives 0.726 / 0.745 / 0.762 —
  never beating ESMFold2 alone. Adding every metric from both methods makes it *worse*
  (hPD-1 0.710 vs 0.757 for ESMFold2's own metrics, p = 0.007): 11 features on 16 positives
  is overfitting, and CV shows it.
- *Structural consensus*: whether the two methods place the binder in the same pose
  (median cross-method RMSD 2.6 Å) or on the same epitope (median Jaccard 0.83) carries
  essentially no signal — AUC 0.46–0.61 across tasks. Two models agreeing is not evidence
  the design binds.

At matched screening budget the ranking barely moves: at k = 25 on hPD-1, ESMFold2 alone
finds 9 binders, fused finds 7, splitting the budget between the two methods finds 8, and an
AND gate finds 5. The AND gate is consistently the worst option — it spends budget on
agreement rather than on candidates.

## Why "no benefit" is not "proven equivalent"

The two methods are only moderately correlated (Spearman ρ = 0.56 human, 0.66 mouse on ipTM),
and they identify *different* binders: at a top-25% cutoff on hPD-1, ESMFold2 flags 9 of the
16 true binders, Protenix-v2 flags 6, and only 5 are flagged by both. So the raw material for
a useful ensemble is there — this dataset just cannot resolve it. With 9–16 positives, the
minimum detectable ΔAUC at 80% power is 0.10–0.21 depending on the task (`Analysis/data/power.csv`).
Every fusion effect measured here is smaller than that.

What would settle it: more positives. Either more designs with BLI data, or a
harder-negative set, or quantitative affinity instead of a binary call at one concentration.

## Sanity checks that passed

- **Binder length is not a confound** — AUC 0.44–0.50 on all three tasks.
- **The natural-ligand control is out of the dataset.** HuPDL1-Fc is labelled `Cross-reactive`
  in the BLI sheet, but it is human PD-L1: 227 aa, a deposited complex. `build_dataset.py`
  drops it from `master_long.csv` and both `merged_*.csv` and parks its row in
  `Analysis/data/positive_control.csv`, so it cannot reach any analysis. Including it had
  inflated the original ESMFold2 figure from 0.758 to 0.77. For the record, from that parked
  row: **Protenix-v2 fails on it outright** (ipTM 0.29 / 0.30, pLDDT 0.55 / 0.63, 50 Å from the
  ESMFold2 pose), while ESMFold2 handles it normally (ipTM 0.86 / 0.84).
- **Independent reimplementation agrees.** `Analysis/plot.py`, patched to read the merged
  tables, reproduces 0.758 / 0.749 / 0.640 / 0.746 and cross-reactive 0.819 / 0.661 exactly.

## Practical recommendation

1. Rank on **ESMFold2 ipTM** (or interface pLDDT — same thing statistically, and it is the
   better of the two on hPD-1).
2. Do **not** spend a Protenix-v2 run to re-rank the same designs; on this evidence it adds
   nothing and on hPD-1 it dilutes.
3. Do **not** use cross-method agreement as a filter.
4. Expect ~2–3× enrichment. Ordering the top 25% of a library recovers a bit over half the
   binders; there is no threshold here that gives a clean pass/fail.
5. AF3 and OpenDDE are still empty. The pipeline handles them the moment their
   `reprediction_metrics.csv` lands — no code change needed.

## Files

| Path | What |
|---|---|
| `Analysis/build_dataset.py` | merges the 4 metric CSVs, recomputes structural features from the 384 complex PDBs |
| `Analysis/predictive_power.py` | AUC / partial AUC / enrichment / DeLong / CV models / power |
| `Analysis/plot_predictive_power.py` | figures 1–6 |
| `Analysis/plot.py` | the original per-method figure, patched for both methods + seaborn 0.11 |
| `Analysis/data/master_long.csv` | tidy: one row per design × species × method (380 rows, 95 designs) |
| `Analysis/data/positive_control.csv` | the removed HuPDL1-Fc rows, parked for reference |
| `Analysis/data/merged_{human,mouse}.csv` | wide per-species tables, `plot.py`-compatible |
| `Analysis/data/auc_table.csv` | every predictor × task with CI, partial AUC, P@K, EF |
| `Analysis/data/head_to_head.csv` | DeLong tests |
| `Analysis/data/budget_matched.csv` | binders found per strategy at fixed budget |
| `Analysis/data/binder_ranks.csv` | every true binder's percentile in each method |
| `Analysis/data/analysis_log.txt` | full console output of the last run |

### Reproduce

```bash
python3 Analysis/build_dataset.py
python3 Analysis/plot_predictive_power.py          # runs the analysis and writes figs 1-6
python3 Analysis/plot.py --human-csv Analysis/data/merged_human.csv \
    --mouse-csv Analysis/data/merged_mouse.csv \
    --classification-csv PD1_duo_DesignSummary.csv --mode both \
    --out Analysis/figures/method_benchmark.png \
    --out-csv Analysis/data/method_benchmark_auc.csv
```

### One thing to fix

`Analysis/plot.py` labels the species axis `PDL1`; the target is PD-1. Two `ax.set_xlabel`
calls, left alone in case it is deliberate shorthand.

---

Rendered as a shareable page: `Analysis/report.html` (self-contained, figures embedded;
regenerate with `python3 Analysis/build_report_html.py`).
