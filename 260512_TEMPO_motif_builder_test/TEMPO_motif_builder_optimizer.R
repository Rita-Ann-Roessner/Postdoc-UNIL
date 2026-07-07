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

objective <- function(mut_weight, decay_factor, vj_prior_strength, len_prior_strength) {
  run_id  <- sprintf("mut%.3f_df%.3f_vj%.1f_len%.1f",
                     mut_weight, decay_factor, vj_prior_strength, len_prior_strength)
  message(sprintf("\n>>> Evaluating: %s", run_id))

  res_list <- list()
  auc01_vals <- sapply(epitopes, function(epitope) {
    mhc     <- sub("_.*", "", epitope)
    peptide <- sub("^[^_]*_", "", epitope)
    out_dir <- file.path(OPTIMIZER_OUTPUT_DIR, run_id, epitope)

    res <- tryCatch(
      run_motif_builder(
        peptide            = peptide,
        mhc                = mhc,
        base_output_dir    = out_dir,
        cdr3_baseline      = cdr3_baseline,
        known_binders_file = file.path(BASE_OUTPUT_DIR, epitope, paste0(epitope, ".csv")),
        validation_file    = file.path(BASE_OUTPUT_DIR, epitope, "validation.csv"),
        mut_weight         = mut_weight,
        decay_factor       = decay_factor,
        vj_baseline_prior  = vj_baseline_prior,
        len_baseline_prior = len_baseline_prior,
        vj_prior_strength  = vj_prior_strength,
        len_prior_strength = len_prior_strength,
        plot_motif         = FALSE,
        plot_final_motif   = TRUE,
        validate_each_step = TRUE
      ),
      error = function(e) {
        message("  ERROR for ", epitope, ": ", e$message)
        NULL
      }
    )

    res_list[[epitope]] <<- res
    if (is.null(res)) NA_real_ else res$final_auc01
  })

  score <- mean(auc01_vals, na.rm = TRUE)
  message(sprintf("    AUC0.1 per epitope: %s  →  mean = %.4f",
                  paste(sprintf("%.4f", auc01_vals), collapse = ", "), score))

  # build long-format rows: one per epitope
  new_rows <- do.call(rbind, lapply(epitopes, function(ep) {
    res_p  <- res_list[[ep]]
    row <- data.frame(
      run_id             = run_id,
      mut_weight         = mut_weight,
      decay_factor       = decay_factor,
      vj_prior_strength  = vj_prior_strength,
      len_prior_strength = len_prior_strength,
      mean_auc01        = score,
      epitope      = ep,
      auc01        = auc01_vals[[ep]],
      stringsAsFactors = FALSE
    )
    # add per-step auc01 as wide columns (auc01_step0, auc01_step1, ...)
    # Pre-allocate all N_ITERATIONS+1 step columns as NA so every epitope's
    # row has the same column set, even if it stopped early (too few top
    # binders) and never reached the later steps.
    for (s in 0:N_ITERATIONS) {
      row[[sprintf("auc01_step%d", s)]] <- NA_real_
    }
    if (!is.null(res_p) && !is.null(res_p$step_counts) &&
        "auc01" %in% colnames(res_p$step_counts)) {
      sc <- res_p$step_counts
      for (i in seq_len(nrow(sc))) {
        col <- sprintf("auc01_step%d", sc$step[i])
        row[[col]] <- sc$auc01[i]
      }
    }
    row
  }))  # end lapply over epitopes

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
    mut_weight         = c(0,    0.3),
    decay_factor       = c(0.05, 0.9),
    vj_prior_strength  = c(0,    100),
    len_prior_strength = c(0,    100)
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
message(sprintf("  MUT_WEIGHT         = %.3f", best["mut_weight"]))
message(sprintf("  DECAY_FACTOR       = %.3f", best["decay_factor"]))
message(sprintf("  VJ_PRIOR_STRENGTH  = %.3f", best["vj_prior_strength"]))
message(sprintf("  LEN_PRIOR_STRENGTH = %.3f", best["len_prior_strength"]))
schedule <- INIT_PERC_RANK * best["decay_factor"]^(0:N_ITERATIONS)
message(sprintf("  threshold schedule: %s",
                paste(sprintf("%.3f", schedule), collapse = " -> ")))
message(sprintf("  Best AUC0.1  = %.4f", bayes_result$Best_Value))
