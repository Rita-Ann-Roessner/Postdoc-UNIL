# =============================================================================
# PAIR_PERCENTILE sweep
#
# Sweeps the selection cut-off applied to an ALREADY FOLDED bait/prey run
# (ESM_motif_builder_pair.R, STEP <- "pair-select" then the cluster fold).
# Nothing is refolded: each cut-off re-selects from the same output_pairs.csv,
# builds its motif and runs TEMPO.
#
# The cut-off is a percentile of the POOLED score distribution over every folded
# pair, both directions together, exactly as FINAL_PERCENTILE cuts the pooled
# alpha+beta distribution of the last step.
#
# Output per cut-off:
#     <BASE>/<epitope>/pairs/p<NN>/
# holding top_pairs.csv, motif_top_binders.csv and validation/. The combined
# table goes to
#     <BASE>/pair_percentile_sweep.csv
# and carries the per-direction diagnostics (additive R^2, bait concordance,
# rho against the dummy ranking, median per-prey sd) alongside the AUCs. Those
# describe the folded matrix rather than the cut-off, so they are computed once
# per epitope and repeat down that epitope's rows.
#
# This leaves the fixed-percentile run that ESM_motif_builder_pair.R writes
# straight into <BASE>/<epitope>/pairs/ untouched, so the two never collide.
#
# Set the config block below.
#
# Usage:  Rscript sweep_pair_percentile.R
# =============================================================================

.sourced_for_benchmark <- TRUE
source("ESM_motif_builder_pair.R")   # helpers + config; main block skipped

# ---- config -----------------------------------------------------------------
PERCENTILES <- c(50, 75, 90, 95, 99)
MOTIF_MODE  <- PAIR_MOTIF_MODE   # "paired" or "chains"; see ESM_motif_builder_pair.R
PLOT_MOTIF  <- FALSE             # MixTCRviz motifs per cut-off, off by default here.
                                 # The fixed-percentile run still writes them.

SUMMARY_CSV <- file.path(BASE_OUTPUT_DIR, "pair_percentile_sweep.csv")

# ---- sweep ------------------------------------------------------------------
rows <- list()

for (label in epitopes) {
  mhc     <- sub("_.*", "", label)
  peptide <- sub("^[^_]*_", "", label)

  scored <- tryCatch(load_pair_scores(label, BASE_OUTPUT_DIR),
                     error = function(e) { warning(sprintf("[%s] %s", label,
                                            conditionMessage(e))); NULL })
  if (is.null(scored)) next

  # matrices + per-prey scores + structure diagnostics, computed once per epitope
  diag <- pair_diagnostics(scored, label, BASE_OUTPUT_DIR)

  validation_file <- file.path(BASE_OUTPUT_DIR, label, "validation.csv")
  if (!(isTRUE(VALIDATION_AVAILABLE) && file.exists(validation_file)))
    message(sprintf("[%s] no validation.csv - motifs will be written, AUC left NA.", label))

  for (p in PERCENTILES) {
    rows[[length(rows) + 1]] <- pair_validate_at(
      scored, percentile = p,
      out_dir = file.path(BASE_OUTPUT_DIR, label, PAIR_DIR, sprintf("p%g", p)),
      label = label, peptide = peptide, mhc = mhc,
      motif_mode = MOTIF_MODE, validation_file = validation_file,
      diag = diag, plot_motif = PLOT_MOTIF)
  }
}

# ---- summary ----------------------------------------------------------------
if (length(rows) == 0) {
  stop("Nothing was scored. Check that ", BASE_OUTPUT_DIR,
       " has a folded <epitope>/", PAIR_DIR, "/output_pairs.csv.")
}
summary_df <- do.call(rbind, rows)
summary_df <- summary_df[order(summary_df$epitope, summary_df$percentile), ]
write.csv(summary_df, SUMMARY_CSV, row.names = FALSE)

message("\n===== PAIR_PERCENTILE sweep =====")
print(summary_df[, c("epitope", "percentile", "threshold", "n_pairs",
                     "n_A_bait", "n_B_bait", "n_alpha", "n_beta",
                     "n_motif", "auc", "auc01")], row.names = FALSE)

message("\n  diagnostics (per fished chain; constant across cut-offs):")
print(unique(summary_df[, c("epitope", "r2_alpha_prey", "concord_alpha",
                            "rho_dummy_alpha", "sd_alpha", "r2_beta_prey",
                            "concord_beta", "rho_dummy_beta", "sd_beta")]),
      row.names = FALSE)

if (any(!is.na(summary_df$auc))) {
  message("\nBest percentile per epitope (by AUC):")
  for (ep in unique(summary_df$epitope)) {
    d <- summary_df[summary_df$epitope == ep & !is.na(summary_df$auc), ]
    if (nrow(d) == 0) next
    b <- d[which.max(d$auc), ]
    message(sprintf("  %-18s p%-4g AUC=%.4f  AUC01=%.4f  (%d alpha, %d beta)",
                    b$epitope, b$percentile, b$auc, b$auc01,
                    b$n_alpha, b$n_beta))
  }
}

message(sprintf("\nWritten to: %s", SUMMARY_CSV))
message(sprintf("Per-cutoff motifs and TEMPO output: %s/<epitope>/%s/p<NN>/",
                BASE_OUTPUT_DIR, PAIR_DIR))
