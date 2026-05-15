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

objective <- function(pssm_weight, decay_factor) {
  run_id  <- sprintf("pw%.2f_df%.3f", pssm_weight, decay_factor)
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
        pssm_weight        = pssm_weight,
        decay_factor       = decay_factor,
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
    data.frame(
      run_id       = run_id,
      pssm_weight  = pssm_weight,
      decay_factor = decay_factor,
      mean_auc01   = score,
      peptide      = p,
      auc01        = auc01_vals[[p]],
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
    pssm_weight  = c(0,    1),
    decay_factor = c(0.05, 0.9)
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
message(sprintf("  PSSM_WEIGHT  = %.3f", best["pssm_weight"]))
message(sprintf("  DECAY_FACTOR = %.3f", best["decay_factor"]))
schedule <- INIT_PERC_RANK * best["decay_factor"]^(0:N_ITERATIONS)
message(sprintf("  threshold schedule: %s",
                paste(sprintf("%.3f", schedule), collapse = " -> ")))
message(sprintf("  Best AUC0.1  = %.4f", bayes_result$Best_Value))
