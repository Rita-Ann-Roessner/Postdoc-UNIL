# =============================================================================
# FINAL_PERCENTILE sweep
#
# Re-runs the "final" selection + TEMPO validation at several percentile cutoffs
# on an ALREADY FOLDED run, so you can see how the final AUC depends on how
# aggressively the last step's top binders are selected.
#
# Nothing is refolded and nothing in the step directories is touched: each
# cutoff gets its own directory
#     <BASE>/<epitope>/validation/p<NN>/
# holding top_tcrs_{alpha,beta}.csv, motif_top_binders.csv, motif_{alpha,beta}/
# and the TEMPO output. A combined table is written to
#     <BASE>/final_percentile_sweep.csv
#
# Usage:  Rscript sweep_final_percentile.R
# =============================================================================

.sourced_for_benchmark <- TRUE
source("ESM_motif_builder.R")   # helpers + config; main block skipped

# ---- config -----------------------------------------------------------------
BASE        <- "TCR_motif_atlas_cdr3_templated_new_installation_reduced_cycles"
PERCENTILES <- c(20, 50, 60, 75, 80, 90, 95)
STEP_INDEX  <- NULL        # NULL = use the highest folded step found per epitope
SUMMARY_CSV <- file.path(BASE, "final_percentile_sweep.csv")

# ---- helpers ----------------------------------------------------------------

# highest step index under <BASE>/<label> that has BOTH chains folded
last_folded_step <- function(label) {
  dirs <- list.files(file.path(BASE, label), pattern = "^step[0-9]+$")
  if (length(dirs) == 0) return(NA_integer_)
  idx <- sort(as.integer(sub("^step", "", dirs)), decreasing = TRUE)
  for (i in idx) {
    d <- file.path(BASE, label, sprintf("step%d", i))
    if (file.exists(file.path(d, "output_alpha.csv")) &&
        file.exists(file.path(d, "output_beta.csv"))) return(i)
  }
  NA_integer_
}

# scored table for one chain of one step dir, or NULL
load_chain <- function(step_dir, chain_letter) {
  chain_name <- if (chain_letter == "A") "alpha" else "beta"
  f <- file.path(step_dir, sprintf("output_%s.csv", chain_name))
  if (!file.exists(f)) return(NULL)
  load_esm_scores(f, meta_csv(step_dir, chain_name))
}

# ---- sweep ------------------------------------------------------------------
rows <- list()

for (label in epitopes) {
  mhc     <- sub("_.*", "", label)
  peptide <- sub("^[^_]*_", "", label)

  step_i <- if (is.null(STEP_INDEX)) last_folded_step(label) else STEP_INDEX
  if (is.na(step_i)) {
    warning(sprintf("[%s] no folded step found under %s - skipped.", label, BASE)); next
  }
  step_dir <- file.path(BASE, label, sprintf("step%d", step_i))

  sA <- load_chain(step_dir, "A")
  sB <- load_chain(step_dir, "B")
  if (is.null(sA) && is.null(sB)) {
    warning(sprintf("[%s] step%d has no usable scores - skipped.", label, step_i)); next
  }

  pooled <- c(if (!is.null(sA)) sA[[SCORE_COL]], if (!is.null(sB)) sB[[SCORE_COL]])
  pooled <- pooled[!is.na(pooled)]

  validation_file <- file.path(BASE, label, "validation.csv")
  have_val <- isTRUE(VALIDATION_AVAILABLE) && file.exists(validation_file)
  if (!have_val)
    message(sprintf("[%s] no validation.csv - motifs will be written, AUC left NA.", label))

  message(sprintf("\n===== %s  (step%d, pooled n=%d) =====", label, step_i, length(pooled)))

  for (pct in PERCENTILES) {
    thr     <- as.numeric(quantile(pooled, pct / 100, na.rm = TRUE))
    out_dir <- file.path(BASE, label, "validation", sprintf("p%g", pct))
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    take <- function(scored, chain_letter) {
      if (is.null(scored)) return(NULL)
      top <- scored[!is.na(scored[[SCORE_COL]]) & scored[[SCORE_COL]] >= thr, ]
      if (nrow(top) == 0) return(NULL)
      chain_name <- if (chain_letter == "A") "alpha" else "beta"
      write.csv(top, file.path(out_dir, sprintf("top_tcrs_%s.csv", chain_name)),
                row.names = FALSE)
      top
    }
    tA <- take(sA, "A"); tB <- take(sB, "B")
    nA <- if (is.null(tA)) 0L else nrow(tA)
    nB <- if (is.null(tB)) 0L else nrow(tB)

    # motif logos, same guard the enrichment steps use
    if (!is.null(tA) && nA >= 10) plot_step_motif(tA, "A", label, step_i, out_dir)
    if (!is.null(tB) && nB >= 10) plot_step_motif(tB, "B", label, step_i, out_dir)

    motif_df <- build_motif_df(tA, tB, paste0(mhc, "_", peptide))
    auc <- auc01 <- NA_real_
    n_motif <- 0L
    if (is.null(motif_df)) {
      warning(sprintf("[%s | p%g] no top binders at threshold %.4f.", label, pct, thr))
    } else {
      n_motif    <- nrow(motif_df)
      motif_file <- file.path(out_dir, "motif_top_binders.csv")
      write.csv(motif_df, motif_file, row.names = FALSE)
      if (have_val) {
        res <- tryCatch(
          validate_motif(motif_file = motif_file, validation_file = validation_file,
                         output_dir = out_dir, peptide = peptide, mhc = mhc),
          error = function(e) { warning(sprintf("[%s | p%g] TEMPO failed: %s",
                                                label, pct, conditionMessage(e))); NULL })
        if (!is.null(res)) { auc <- res$auc; auc01 <- res$auc01 }
      }
    }

    message(sprintf("  p%-4g thr=%.4f  alpha=%-5d beta=%-5d motif=%-5d AUC=%s",
                    pct, thr, nA, nB, n_motif,
                    if (is.na(auc)) "NA" else sprintf("%.4f", auc)))

    rows[[length(rows) + 1]] <- data.frame(
      epitope = label, step = step_i, percentile = pct, threshold = round(thr, 4),
      n_alpha = nA, n_beta = nB, n_motif = n_motif,
      auc = auc, auc01 = auc01, stringsAsFactors = FALSE)
  }
}

# ---- summary ----------------------------------------------------------------
if (length(rows) == 0) {
  stop("Nothing was scored. Check that ", BASE, " has folded step directories.")
}
summary_df <- do.call(rbind, rows)
summary_df <- summary_df[order(summary_df$epitope, summary_df$percentile), ]
write.csv(summary_df, SUMMARY_CSV, row.names = FALSE)

message("\n===== FINAL_PERCENTILE sweep =====")
print(summary_df, row.names = FALSE)
message(sprintf("\nWritten to: %s", SUMMARY_CSV))
message(sprintf("Per-cutoff motifs and TEMPO output: %s/<epitope>/validation/p<NN>/", BASE))
