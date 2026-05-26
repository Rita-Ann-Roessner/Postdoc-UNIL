# TEMPO Motif Builder

An iterative pipeline for generating synthetic TCR pools enriched for epitope-specific binding, using [TEMPOtrain](https://github.com/GfellerLab/TEMPO) as a scoring backbone and [MixTCRviz](https://github.com/GfellerLab/MixTCRviz) for motif visualization.

## Overview

The pipeline starts from a fully random TCR pool and iteratively enriches it toward TCRs predicted to bind a given peptide–MHC complex. At each step, the top-scoring TCRs are used to bias the V/J gene usage, CDR3 length distribution, and CDR3 amino acid composition (via a PSSM) of the next generated pool. Performance is evaluated by AUC0.1 on a labelled validation set.

### Pipeline steps

```
Step 0   Random TCR generation (flat V/J baseline) → TEMPO scoring → top binders
Step 1–N Enrich from top binders (V/J + CDR3 length + PSSM) → generate → score → top binders
Final    Validate last model on labelled validation.csv → AUC / AUC0.1
```

### Three enrichment signals

| Signal | Description |
|--------|-------------|
| **V/J joint distribution** | Sample (V,J) pairs proportionally to their frequency in top binders |
| **CDR3 length** | Re-weight length sampling toward lengths enriched in top binders (marginal or conditioned on V/J) |
| **PSSM** | Blend per-position amino acid probabilities between V/J baseline and a motif built from top binders |

## Files

| File | Description |
|------|-------------|
| `TEMPO_motif_builder_test.R` | Main pipeline — functions + config + run section |
| `TEMPO_motif_builder_optimizer.R` | Bayesian optimisation of `pssm_weight` and `decay_factor` |
| `TEMPO_grid_search.R` | Grid search over `n_pairs` × `n_cdr3_multi` |
| `build_reference_motifs.R` | Builds reference MixTCRviz motifs from known binders |
| `Monitor.ipynb` | Result visualisation (AUC curves, grid search plots) |
| `test.R` | TEMPO trained on known binders (TCR motif atlas data) → upper ceiling |

## Configuration

Key parameters in `TEMPO_motif_builder_test.R`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `N_ITERATIONS` | 3 | Number of enrichment steps (steps 1–N) |
| `N_PAIRS` | 400 | V/J pairs sampled per chain per step |
| `N_CDR3_MULTI` | 5 | CDR3 sequences sampled per V/J pair |
| `INIT_PERC_RANK` | 10 | Percentile rank threshold at step 0 |
| `DECAY_FACTOR` | 0.135 | Multiplicative threshold decay per step |
| `PSSM_WEIGHT` | 0.643 | Blend weight between baseline and PSSM (0 = pure baseline, 1 = pure PSSM) |
| `MIN_TCRS_PSSM` | 30 | Minimum top binders required before building a PSSM |
| `LEN_DIST_COND_VJ` | FALSE | If TRUE, sample CDR3 length conditioned on V/J pair |

The threshold at step *k* follows a geometric decay schedule:  
`threshold(k) = INIT_PERC_RANK × DECAY_FACTOR^k`

## Input data

Each epitope lives in its own subdirectory, e.g. `ELAGIGILTV/`:

```
{peptide}/
  A0201_{peptide}.csv   # known binders (used to build TEMPO predictor + reference motif)
  validation.csv        # labelled TCRs for AUC evaluation (Label column required)
```

## Output structure

```
predictor_{peptide}.rds          # cached TEMPO predictor (built once, reused)
TEMPO_step0_{peptide}/           # step 0: random pool
  model.csv                      # scored TCR pool
  top_tcrs.csv                   # TCRs passing the step threshold
  TEMPO_scoring/                 # raw TEMPO output
TEMPO_step{1..N}_{peptide}/      # enrichment steps
  model.csv / top_tcrs.csv / TEMPO_scoring/
  motif_vs_iter/                 # MixTCRviz motif vs all iterated TCRs (final step only)
  motif_vs_baseline/             # MixTCRviz motif vs standard baseline (final step only)
TEMPO_validation_{peptide}/      # final AUC validation output
{peptide}/motif_plain.png        # reference motif from known binders
optimizer_runs/                  # Bayesian optimisation results
grid_search_runs/                # Grid search results
  grid_results.csv               # AUC / AUC0.1 per peptide × parameter combination
```

## Epitopes

The pipeline has been tested on four HLA-A*02:01-restricted epitopes:

| Epitope | Source |
|---------|--------|
| ELAGIGILTV | MART-1 |
| GILGFVFTL | Influenza M1 |
| GLCTLVAML | EBV BMLF1 |
| LLWNGPMAV | HIV Gag |

## Dependencies

```r
install.packages(c("TEMPOtrain", "MixTCRviz", "dplyr", "tidyr", "pROC",
                   "rBayesianOptimization", "pdftools"))
```
