# =============================================================================
# Benchmark generator — VJ_PRIOR_STRENGTH sweep (coordinate-descent: VJ first,
# LEN fixed). Generates the step-<STEP> ESMFold input batches for every grid
# point, each branching from a shared, already-scored previous step copied from
# TCR_motif_atlas. Only the enrichment steps (>=1) depend on the knobs, so step 0
# is shared across all grid points.
#
# Round-trip workflow (STEP given on the command line, default 1):
#   Rscript benchmark_generate.R 1     # seed step0, make step1 batches
#   ./benchmark_submit.sh 1  ; (wait) ; ./benchmark_gather.sh 1     # fold on cluster
#   Rscript benchmark_generate.R 2     # make step2 batches from folded step1
#   ./benchmark_submit.sh 2  ; (wait) ; ./benchmark_gather.sh 2
#   Rscript benchmark_evaluate.R       # AUC + false-positive-gene trajectory
# =============================================================================

.sourced_for_benchmark <- TRUE
source("ESM_motif_builder.R")   # functions + setup (priors, gene tables); main loop skipped

# ---- Benchmark configuration ------------------------------------------------
BENCH_DIR     <- "benchmark_MUT_v2"
# Sweep VJ, LEN, or both: a vector sweeps that knob, a scalar fixes it. When both
# are vectors the double loop below generates the full VJ x LEN grid.
GRID_VJ       <- 60                     # VJ_PRIOR_STRENGTH (fixed at the VJ-sweep winner)
GRID_LEN      <- 20          # LEN_PRIOR_STRENGTH sweep
BENCH_N_PAIRS <- 400                    # smaller than the full-run N_PAIRS to cut folds
SRC_DIR       <- BASE_OUTPUT_DIR        # "TCR_motif_atlas" — source of the shared step0

STEP <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)[1]))
if (is.na(STEP)) STEP <- 2

gp_name <- function(vj, len) sprintf("vj%03d_len%03d", vj, len)

for (vj in GRID_VJ) for (len in GRID_LEN) {
  gp     <- gp_name(vj, len)
  gp_dir <- file.path(BENCH_DIR, gp)

  for (ep in epitopes) {
    mhc        <- sub("_.*", "", ep)
    peptide    <- sub("^[^_]*_", "", ep)
    mhc_allele <- if (grepl("^HLA_", mhc)) mhc else paste0("HLA_", mhc)

    # On the first generation (step 1), seed each grid tree with the SHARED,
    # knob-independent step0 (model + scores) and the validation set.
    if (STEP == 1) {
      s0_dst <- file.path(gp_dir, ep, "step0")
      dir.create(s0_dst, recursive = TRUE, showWarnings = FALSE)
      for (f in c("model_alpha.csv", "model_beta.csv",
                  "output_alpha.csv", "output_beta.csv")) {
        src <- file.path(SRC_DIR, ep, "step0", f)
        if (!file.exists(src)) stop("missing shared step0 file: ", src)
        file.copy(src, file.path(s0_dst, f), overwrite = TRUE)
      }
      file.copy(file.path(SRC_DIR, ep, "validation.csv"),
                file.path(gp_dir, ep, "validation.csv"), overwrite = TRUE)
    }

    message(sprintf("\n[grid %s | %s] generating step %d (VJ=%d, LEN=%d)",
                    gp, ep, STEP, vj, len))
    run_enrich_step(
      step               = STEP,
      peptide            = peptide,
      mhc_allele         = mhc_allele,
      label              = ep,
      cdr3_baseline      = cdr3_baseline,
      base_output_dir    = gp_dir,
      n_pairs            = BENCH_N_PAIRS,
      n_cdr3_multi       = N_CDR3_MULTI,
      min_tcrs_pssm      = MIN_TCRS_PSSM,
      vj_baseline_prior  = vj_baseline_prior,
      vj_prior_strength  = vj,            # <-- the swept knob
      len_baseline_prior = len_baseline_prior,
      len_prior_strength = len,
      esm_thresholds     = ESM_THRESHOLDS,
      plot_each_step     = FALSE,
      validate_each_step = FALSE,
      validation_file    = file.path(gp_dir, ep, "validation.csv"),
      mhc                = mhc
    )
  }
}

message(sprintf("\nDone. Fold on the cluster, e.g.:\n  ./benchmark_submit.sh %d\n  (wait)\n  ./benchmark_gather.sh %d",
                STEP, STEP))
