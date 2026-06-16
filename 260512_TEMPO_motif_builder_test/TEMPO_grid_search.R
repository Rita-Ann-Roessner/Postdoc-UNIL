# =============================================================================
# Grid search over VJ_PRIOR_STRENGTH
# Sources TEMPO_motif_builder_test.R (functions + config only, run section
# is skipped via the .sourced_by_optimizer guard).
# Saves AUC / AUC0.1 per peptide × parameter combination to a CSV.
# =============================================================================

.sourced_by_optimizer <- TRUE
source("TEMPO_motif_builder_test.R", local = FALSE)

GRID_OUTPUT_DIR <- "grid_search_runs_prior"
dir.create(GRID_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Grid definition --------------------------------------------------------

grid <- expand.grid(
  vj_prior_strength = c(0, 5, 20, 50, 100),
  stringsAsFactors = FALSE
)

# ---- Run --------------------------------------------------------------------

results <- data.frame()

for (i in seq_len(nrow(grid))) {

  vps <- grid$vj_prior_strength[i]
  run_id <- sprintf("vjprior%.1f", vps)
  message(sprintf("\n===== Grid run %d / %d : %s =====", i, nrow(grid), run_id))

  for (epitope in epitopes) {
    mhc     <- sub("_.*", "", epitope)
    peptide <- sub("^[^_]*_", "", epitope)
    out_dir <- file.path(GRID_OUTPUT_DIR, run_id, epitope)
    message(sprintf("  -> %s", epitope))

    res <- tryCatch(
      run_motif_builder(
        peptide            = peptide,
        mhc                = mhc,
        base_output_dir    = out_dir,
        cdr3_baseline      = cdr3_baseline,
        known_binders_file = file.path(BASE_OUTPUT_DIR, epitope, paste0(epitope, ".csv")),
        validation_file    = file.path(BASE_OUTPUT_DIR, epitope, "validation.csv"),
        vj_baseline_prior  = vj_baseline_prior,
        vj_prior_strength  = vps,
        plot_motif         = FALSE,
        plot_final_motif   = FALSE,
        validate_each_step = FALSE
      ),
      error = function(e) {
        message("  ERROR: ", e$message)
        NULL
      }
    )

    row <- data.frame(
      run_id            = run_id,
      vj_prior_strength = vps,
      epitope           = epitope,
      auc               = if (is.null(res)) NA_real_ else res$final_auc,
      auc01             = if (is.null(res)) NA_real_ else res$final_auc01,
      stringsAsFactors = FALSE
    )
    results <- rbind(results, row)

    # write incrementally so results are not lost on early exit
    write.csv(results,
              file.path(GRID_OUTPUT_DIR, "grid_results.csv"),
              row.names = FALSE)
  }
}

# ---- Summary ----------------------------------------------------------------

message("\n===== Grid search complete =====")
message(sprintf("Results saved to: %s", file.path(GRID_OUTPUT_DIR, "grid_results.csv")))
print(results)
