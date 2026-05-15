library(rBayesianOptimization)

# Source all functions and config from the main script.
# The ### ---- Run ---- section at the bottom will also execute (building
# baseline, dico, etc.) but NOT run run_motif_builder — that loop is in the
# run section which we skip by sourcing with local=TRUE and then detaching.
# Simplest approach: source the whole script but wrap the run block in an
# if(FALSE) guard there, OR just source functions only by stopping before the
# run section. We use sys.source with a max-line trick below.

.sourced_by_optimizer <- TRUE
source("TEMPO_motif_builder_test.R", local = FALSE)

# Working directory should be:
# /Users/roessner/Documents/PostDoc/Data/260512_TEMPO_motif_builder_test

OPTIMIZER_OUTPUT_DIR <- "optimizer_runs"
dir.create(OPTIMIZER_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Objective function
# Maximise mean AUC0.1 across all three epitopes.
# N_ENRICHMENT is continuous here — rounded to integer inside.
# plot_motif=FALSE suppresses MixTCRviz PDF generation during search.
# =============================================================================

eval_log <- data.frame()

objective <- function(n_enrichment, pssm_weight, perc_rank_mid) {
  n_enrichment <- as.integer(round(n_enrichment))

  run_id  <- sprintf("ne%02d_pw%.2f_prs%.2f", n_enrichment, pssm_weight, perc_rank_mid)
  message(sprintf("\n>>> Evaluating: %s", run_id))

  res_list <- list()
  auc01_vals <- sapply(names(dico), function(peptide) {
    mhc     <- dico[[peptide]]
    out_dir <- file.path(OPTIMIZER_OUTPUT_DIR, run_id, peptide)

    res <- tryCatch(
      run_motif_builder(
        peptide            = peptide,
        mhc                = mhc,
        base_output_dir    = out_dir,
        cdr3_baseline      = cdr3_baseline,
        known_binders_file = file.path(peptide, paste0("A0201_", peptide, ".csv")),
        validation_file    = file.path(peptide, "validation.csv"),
        n_enrichment       = n_enrichment,
        pssm_weight        = pssm_weight,
        perc_rank_mid    = perc_rank_mid,
        plot_motif         = FALSE
      ),
      error = function(e) {
        message("  ERROR for ", peptide, ": ", e$message)
        NULL
      }
    )

    res_list[[peptide]] <<- res
    if (is.null(res)) NA_real_ else res$final_auc01
  })

  score <- mean(auc01_vals, na.rm = TRUE)
  message(sprintf("    AUC0.1 per epitope: %s  →  mean = %.4f",
                  paste(sprintf("%.4f", auc01_vals), collapse = ", "), score))

  # build long-format rows: one per peptide × step
  new_rows <- do.call(rbind, lapply(names(dico), function(p) {
    res_p  <- res_list[[p]]
    counts <- if (!is.null(res_p) && !is.null(res_p$step_counts))
                res_p$step_counts
              else
                data.frame(step = NA, n_total = NA, n_top = NA)
    data.frame(
      run_id          = run_id,
      n_enrichment    = n_enrichment,
      pssm_weight     = pssm_weight,
      perc_rank_mid = perc_rank_mid,
      mean_auc01      = score,
      peptide         = p,
      auc01           = auc01_vals[[p]],
      step            = counts$step,
      n_total         = counts$n_total,
      n_top           = counts$n_top,
      stringsAsFactors = FALSE
    )
  }))

  eval_log <<- rbind(eval_log, new_rows)
  write.csv(eval_log,
            file.path(OPTIMIZER_OUTPUT_DIR, "bayes_history.csv"),
            row.names = FALSE)

  list(Score = score, Pred = 0)
}

# =============================================================================
# Run Bayesian optimisation
# init_points: random exploration before fitting the GP surrogate
# n_iter:      Bayesian iterations after init
# =============================================================================

set.seed(42)

bayes_result <- BayesianOptimization(
  FUN         = objective,
  bounds      = list(
    n_enrichment    = c(5L,  20L),
    pssm_weight     = c(0,   1),
    perc_rank_mid = c(1,   10)
  ),
  init_points = 8,
  n_iter      = 20,
  acq         = "ucb",
  kappa       = 2.576,
  verbose     = TRUE
)

# =============================================================================
# Save and display results
# =============================================================================

# eval_log already written incrementally during the run

best <- bayes_result$Best_Par
message("\n===== Best parameters =====")
message(sprintf("  N_ENRICHMENT    = %d",   as.integer(round(best["n_enrichment"]))))
message(sprintf("  PSSM_WEIGHT     = %.3f", best["pssm_weight"]))
message(sprintf("  PERC_RANK_MID = %.3f", best["perc_rank_mid"]))
message(sprintf("  Best AUC0.1     = %.4f", bayes_result$Best_Value))

schedule <- exp(seq(log(best["perc_rank_mid"]), log(0.05), length.out = N_ITERATIONS))
message(sprintf("  perc_rank schedule: 10 → %s",
                paste(sprintf("%.3f", schedule), collapse = " → ")))
