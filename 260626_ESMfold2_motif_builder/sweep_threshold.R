# =============================================================================
# Absolute-threshold sweep
#
# For an ALREADY FOLDED step, re-runs the top-binder selection + motif build +
# TEMPO validation at several FIXED iptm_pair_mean cut-offs, so you can see how
# the motif quality depends on where the threshold sits. Use it to pick
# ESM_THRESHOLDS[step + 1] without folding anything.
#
# Nothing is refolded and nothing in the step directories is touched: each
# cut-off gets its own directory
#     <BASE>/<epitope>/validation/step<N>_t<thr>/
# holding top_tcrs_{alpha,beta}.csv, motif_top_binders.csv, motif_{alpha,beta}/
# and the TEMPO output. A combined table is written to
#     <BASE>/threshold_sweep_step<N>.csv
#
# The threshold is applied to each chain separately, exactly as
# enrich_one_chain does, so the per-chain pass RATES will differ. Both the
# counts and the rates are logged.
#
# Set BASE / STEP_INDEX / THRESHOLDS in the config block below.
#
# Usage:  Rscript sweep_threshold.R
# =============================================================================

.sourced_for_benchmark <- TRUE
source("ESM_motif_builder.R")   # helpers + config; main block skipped

# ---- config -----------------------------------------------------------------
BASE       <- "TCR_motif_atlas_cdr3_templated_new_installation_reduced_cycles"
STEP_INDEX <- 1L          # which folded step to select from; NULL = highest folded
THRESHOLDS <- c(0.55, 0.56, 0.60, 0.65, 0.70, 0.75, 0.80)

SUMMARY_CSV <- file.path(BASE, sprintf("threshold_sweep_step%s.csv",
                                       if (is.null(STEP_INDEX)) "last" else STEP_INDEX))

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

n_scored <- function(scored) {
  if (is.null(scored)) 0L else sum(!is.na(scored[[SCORE_COL]]))
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
  na_tot <- n_scored(sA); nb_tot <- n_scored(sB)

  validation_file <- file.path(BASE, label, "validation.csv")
  have_val <- isTRUE(VALIDATION_AVAILABLE) && file.exists(validation_file)
  if (!have_val)
    message(sprintf("[%s] no validation.csv - motifs will be written, AUC left NA.", label))

  message(sprintf("\n===== %s  (step%d, scored: alpha=%d beta=%d) =====",
                  label, step_i, na_tot, nb_tot))

  for (thr in THRESHOLDS) {
    out_dir <- file.path(BASE, label, "validation",
                         sprintf("step%d_t%.3f", step_i, thr))
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    take <- function(scored, chain_letter) {
      if (is.null(scored)) return(NULL)
      top <- scored[!is.na(scored[[SCORE_COL]]) & scored[[SCORE_COL]] >= thr, ]
      if (nrow(top) == 0) return(NULL)
      chain_name <- if (chain_letter == "A") "alpha" else "beta"
      write.csv(top, file.path(out_dir, sprintf("top_tcrs_%s.csv", chain_name)),
                row.names = FALSE)
      write_fold_input_csv(top, file.path(out_dir,
                           sprintf("top_tcrs_%s_seqs.csv", chain_name)))
      top
    }
    tA <- take(sA, "A"); tB <- take(sB, "B")
    nA <- if (is.null(tA)) 0L else nrow(tA)
    nB <- if (is.null(tB)) 0L else nrow(tB)
    pA <- if (na_tot > 0) 100 * nA / na_tot else NA_real_
    pB <- if (nb_tot > 0) 100 * nB / nb_tot else NA_real_

    # motif logos, same guard the enrichment steps use
    if (!is.null(tA) && nA >= 10) plot_step_motif(tA, "A", label, step_i, out_dir)
    if (!is.null(tB) && nB >= 10) plot_step_motif(tB, "B", label, step_i, out_dir)

    motif_df <- build_motif_df(tA, tB, paste0(mhc, "_", peptide))
    auc <- auc01 <- NA_real_
    n_motif <- 0L
    if (is.null(motif_df)) {
      warning(sprintf("[%s | t%.3f] no top binders at this threshold.", label, thr))
    } else {
      n_motif    <- nrow(motif_df)
      motif_file <- file.path(out_dir, "motif_top_binders.csv")
      write.csv(motif_df, motif_file, row.names = FALSE)
      if (have_val) {
        res <- tryCatch(
          validate_motif(motif_file = motif_file, validation_file = validation_file,
                         output_dir = out_dir, peptide = peptide, mhc = mhc),
          error = function(e) { warning(sprintf("[%s | t%.3f] TEMPO failed: %s",
                                                label, thr, conditionMessage(e))); NULL })
        if (!is.null(res)) { auc <- res$auc; auc01 <- res$auc01 }
      }
    }

    # flag cells too thin to build a PSSM from, whatever their AUC says
    thin <- (nA > 0 && nA < MIN_TCRS_PSSM) || (nB > 0 && nB < MIN_TCRS_PSSM)
    message(sprintf("  t=%.3f  alpha=%-5d (%4.1f%%)  beta=%-5d (%4.1f%%)  motif=%-5d AUC=%s%s",
                    thr, nA, pA, nB, pB, n_motif,
                    if (is.na(auc)) "NA" else sprintf("%.4f", auc),
                    if (thin) "   [< MIN_TCRS_PSSM]" else ""))

    rows[[length(rows) + 1]] <- data.frame(
      epitope = label, step = step_i, threshold = thr,
      n_alpha = nA, n_scored_alpha = na_tot, pct_alpha = round(pA, 1),
      n_beta  = nB, n_scored_beta  = nb_tot, pct_beta  = round(pB, 1),
      n_motif = n_motif, thin_pssm = thin,
      auc = auc, auc01 = auc01, stringsAsFactors = FALSE)
  }
}

# ---- summary ----------------------------------------------------------------
if (length(rows) == 0) {
  stop("Nothing was scored. Check that ", BASE, " has folded step directories.")
}
summary_df <- do.call(rbind, rows)
summary_df <- summary_df[order(summary_df$epitope, summary_df$threshold), ]
write.csv(summary_df, SUMMARY_CSV, row.names = FALSE)

message("\n===== absolute-threshold sweep =====")
print(summary_df, row.names = FALSE)

# best threshold per epitope, by AUC
if (any(!is.na(summary_df$auc))) {
  message("\nBest AUC per epitope:")
  for (ep in unique(summary_df$epitope)) {
    d <- summary_df[summary_df$epitope == ep & !is.na(summary_df$auc), ]
    if (nrow(d) == 0) next
    b <- d[which.max(d$auc), ]
    message(sprintf("  %-18s t=%.3f  AUC=%.4f  (alpha=%d, beta=%d)",
                    b$epitope, b$threshold, b$auc, b$n_alpha, b$n_beta))
  }
}

message(sprintf("\nWritten to: %s", SUMMARY_CSV))
message(sprintf("Per-threshold motifs and TEMPO output: %s/<epitope>/validation/step%s_t<thr>/",
                BASE, if (is.null(STEP_INDEX)) "<N>" else STEP_INDEX))
