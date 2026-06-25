# =============================================================================
# Grid search over PSSM_WEIGHT (w_max)
# Sources TEMPO_motif_builder_test.R (functions + config only, run section
# is skipped via the .sourced_by_optimizer guard).
# Saves AUC / AUC0.1 per peptide × parameter combination to a CSV.
#
# NOTE: the objective is stochastic (random TCR sampling), so each grid point
# is evaluated over N_REPLICATES seeds. This reveals the noise band rather than
# reading too much into a single run.
# =============================================================================

.sourced_by_optimizer <- TRUE
source("TEMPO_motif_builder_test.R", local = FALSE)

GRID_OUTPUT_DIR <- "grid_search_runs_pssm"
dir.create(GRID_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

N_REPLICATES <- 3   # seeds per grid point (noisy objective)

# ---- Grid definition --------------------------------------------------------

grid <- expand.grid(
  pssm_weight = seq(0.3, 1.0, by = 0.1),
  stringsAsFactors = FALSE
)

# ---- Run --------------------------------------------------------------------

results <- data.frame()

for (i in seq_len(nrow(grid))) {

  pw <- grid$pssm_weight[i]
  message(sprintf("\n===== Grid run %d / %d : pssm_weight = %.2f =====",
                  i, nrow(grid), pw))

  for (rep in seq_len(N_REPLICATES)) {
    set.seed(rep)
    run_id <- sprintf("pssm%.2f_rep%d", pw, rep)
    message(sprintf("  -- replicate %d / %d", rep, N_REPLICATES))

    for (epitope in epitopes) {
      mhc     <- sub("_.*", "", epitope)
      peptide <- sub("^[^_]*_", "", epitope)
      out_dir <- file.path(GRID_OUTPUT_DIR, run_id, epitope)
      message(sprintf("    -> %s", epitope))

      res <- tryCatch(
        run_motif_builder(
          peptide            = peptide,
          mhc                = mhc,
          base_output_dir    = out_dir,
          cdr3_baseline      = cdr3_baseline,
          known_binders_file = file.path(BASE_OUTPUT_DIR, epitope, paste0(epitope, ".csv")),
          validation_file    = file.path(BASE_OUTPUT_DIR, epitope, "validation.csv"),
          pssm_weight        = pw,
          vj_baseline_prior  = vj_baseline_prior,
          plot_motif         = FALSE,
          plot_final_motif   = FALSE,
          validate_each_step = FALSE
        ),
        error = function(e) {
          message("    ERROR: ", e$message)
          NULL
        }
      )

      row <- data.frame(
        run_id      = run_id,
        pssm_weight = pw,
        replicate   = rep,
        epitope     = epitope,
        auc         = if (is.null(res)) NA_real_ else res$final_auc,
        auc01       = if (is.null(res)) NA_real_ else res$final_auc01,
        stringsAsFactors = FALSE
      )
      results <- rbind(results, row)

      # write incrementally so results are not lost on early exit
      write.csv(results,
                file.path(GRID_OUTPUT_DIR, "grid_results.csv"),
                row.names = FALSE)
    }
  }
}

# ---- Summary ----------------------------------------------------------------

message("\n===== Grid search complete =====")
message(sprintf("Results saved to: %s", file.path(GRID_OUTPUT_DIR, "grid_results.csv")))

# Mean AUC0.1 per pssm_weight, averaged across epitopes and replicates
agg <- aggregate(auc01 ~ pssm_weight, data = results, FUN = mean, na.rm = TRUE)
message("\nMean AUC0.1 by pssm_weight (across epitopes × replicates):")
print(agg)
print(results)
