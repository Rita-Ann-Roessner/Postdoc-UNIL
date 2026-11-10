# =============================================================================
# CHAIN_WEIGHT sweep
#
# For an ALREADY FOLDED step, selects the top binders once and then runs the
# TEMPO validation three times, with chain_weight = "equal", "sqrt" and
# "linear", so you can see how the alpha/beta combination rule affects the AUC.
#
# Nothing is refolded and nothing in the step directories is touched. The motif
# is identical across the three methods (it depends only on the selection), so
# it is built once per epitope and reused. Each method gets its own directory
#     <BASE>/<epitope>/validation/step<N>_<method>/
# holding validation_pred.csv and the TEMPO output. A combined table goes to
#     <BASE>/chain_weight_sweep_step<N>.csv
#
# Set the config block below.
#
# Usage:  Rscript sweep_chain_weight.R
# =============================================================================

.sourced_for_benchmark <- TRUE
source("ESM_motif_builder.R")   # helpers + config; main block skipped

# ---- config -----------------------------------------------------------------
BASE       <- "TCR_motif_atlas_cdr3_templated_new_installation_reduced_cycles"
STEP_INDEX <- 2L               # which folded step to select the top binders from
METHODS    <- c("equal", "sqrt", "linear")

# How to pick the top binders. Set SELECT_PERCENTILE to a number to use a
# percentile of the POOLED (alpha+beta) score distribution; leave it NULL to use
# the fixed SELECT_THRESHOLD instead.
SELECT_PERCENTILE <- NULL
SELECT_THRESHOLD  <- 0.67

SUMMARY_CSV <- file.path(BASE, sprintf("chain_weight_sweep_step%d.csv", STEP_INDEX))

# ---- helpers ----------------------------------------------------------------

load_chain <- function(step_dir, chain_letter) {
  chain_name <- if (chain_letter == "A") "alpha" else "beta"
  f <- file.path(step_dir, sprintf("output_%s.csv", chain_name))
  if (!file.exists(f)) return(NULL)
  load_esm_scores(f, meta_csv(step_dir, chain_name))
}

n_scored <- function(scored) if (is.null(scored)) 0L else sum(!is.na(scored[[SCORE_COL]]))

# ---- sweep ------------------------------------------------------------------
rows <- list()

for (label in epitopes) {
  mhc     <- sub("_.*", "", label)
  peptide <- sub("^[^_]*_", "", label)

  step_dir <- file.path(BASE, label, sprintf("step%d", STEP_INDEX))
  sA <- load_chain(step_dir, "A")
  sB <- load_chain(step_dir, "B")
  if (is.null(sA) && is.null(sB)) {
    warning(sprintf("[%s] step%d has no usable scores - skipped.", label, STEP_INDEX)); next
  }
  na_tot <- n_scored(sA); nb_tot <- n_scored(sB)

  validation_file <- file.path(BASE, label, "validation.csv")
  if (!(isTRUE(VALIDATION_AVAILABLE) && file.exists(validation_file))) {
    warning(sprintf("[%s] no validation.csv - skipped.", label)); next
  }

  # ---- one selection, shared by all three methods ----
  if (!is.null(SELECT_PERCENTILE)) {
    pooled <- c(if (!is.null(sA)) sA[[SCORE_COL]], if (!is.null(sB)) sB[[SCORE_COL]])
    pooled <- pooled[!is.na(pooled)]
    thr    <- as.numeric(quantile(pooled, SELECT_PERCENTILE / 100, na.rm = TRUE))
    sel_lab <- sprintf("p%g (thr %.4f)", SELECT_PERCENTILE, thr)
  } else {
    thr     <- SELECT_THRESHOLD
    sel_lab <- sprintf("thr %.4f", thr)
  }

  take <- function(scored) {
    if (is.null(scored)) return(NULL)
    top <- scored[!is.na(scored[[SCORE_COL]]) & scored[[SCORE_COL]] >= thr, ]
    if (nrow(top) == 0) NULL else top
  }
  tA <- take(sA); tB <- take(sB)
  nA <- if (is.null(tA)) 0L else nrow(tA)
  nB <- if (is.null(tB)) 0L else nrow(tB)

  motif_df <- build_motif_df(tA, tB, paste0(mhc, "_", peptide))
  if (is.null(motif_df)) {
    warning(sprintf("[%s] no top binders at %s - skipped.", label, sel_lab)); next
  }

  base_dir <- file.path(BASE, label, "validation")
  dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
  motif_file <- file.path(base_dir, sprintf("step%d_chainweight_motif.csv", STEP_INDEX))
  write.csv(motif_df, motif_file, row.names = FALSE)

  # motif-level chain counts, i.e. what the weights are actually derived from
  mA <- sum(!is.na(motif_df$cdr3_TRA)); mB <- sum(!is.na(motif_df$cdr3_TRB))

  message(sprintf("\n===== %s  (step%d, %s) =====", label, STEP_INDEX, sel_lab))
  message(sprintf("  selected alpha=%d/%d (%.1f%%)  beta=%d/%d (%.1f%%)  motif rows=%d",
                  nA, na_tot, if (na_tot) 100 * nA / na_tot else NA_real_,
                  nB, nb_tot, if (nb_tot) 100 * nB / nb_tot else NA_real_,
                  nrow(motif_df)))

  for (method in METHODS) {
    out_dir <- file.path(base_dir, sprintf("step%d_%s", STEP_INDEX, method))
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    wfun <- switch(method, equal = function(x) 1, sqrt = sqrt, linear = identity)
    wA <- wfun(mA); wB <- wfun(mB)

    res <- tryCatch(
      validate_motif(motif_file = motif_file, validation_file = validation_file,
                     output_dir = out_dir, peptide = peptide, mhc = mhc,
                     chain_weight = method),
      error = function(e) { warning(sprintf("[%s | %s] TEMPO failed: %s",
                                            label, method, conditionMessage(e))); NULL })
    auc <- if (is.null(res)) NA_real_ else res$auc
    auc01 <- if (is.null(res)) NA_real_ else res$auc01

    message(sprintf("  %-7s wA=%.2f wB=%.2f  AUC=%s  AUC01=%s",
                    method, wA, wB,
                    if (is.na(auc))   "NA" else sprintf("%.4f", auc),
                    if (is.na(auc01)) "NA" else sprintf("%.4f", auc01)))

    rows[[length(rows) + 1]] <- data.frame(
      epitope = label, step = STEP_INDEX, method = method,
      threshold = round(thr, 4),
      n_alpha = nA, n_scored_alpha = na_tot,
      n_beta  = nB, n_scored_beta  = nb_tot,
      n_motif = nrow(motif_df), motif_alpha = mA, motif_beta = mB,
      w_alpha = wA, w_beta = wB,
      auc = auc, auc01 = auc01, stringsAsFactors = FALSE)
  }
}

# ---- summary ----------------------------------------------------------------
if (length(rows) == 0) {
  stop("Nothing was scored. Check that ", BASE, " has a folded step", STEP_INDEX, ".")
}
summary_df <- do.call(rbind, rows)
summary_df <- summary_df[order(summary_df$epitope, summary_df$method), ]
write.csv(summary_df, SUMMARY_CSV, row.names = FALSE)

message("\n===== CHAIN_WEIGHT sweep =====")
print(summary_df[, c("epitope", "method", "threshold", "n_alpha", "n_beta",
                     "w_alpha", "w_beta", "auc", "auc01")], row.names = FALSE)

if (any(!is.na(summary_df$auc))) {
  message("\nBest method per epitope (by AUC):")
  for (ep in unique(summary_df$epitope)) {
    d <- summary_df[summary_df$epitope == ep & !is.na(summary_df$auc), ]
    if (nrow(d) == 0) next
    b <- d[which.max(d$auc), ]
    message(sprintf("  %-18s %-7s AUC=%.4f  AUC01=%.4f", b$epitope, b$method, b$auc, b$auc01))
  }
}

message(sprintf("\nWritten to: %s", SUMMARY_CSV))
message(sprintf("Per-method TEMPO output: %s/<epitope>/validation/step%d_<method>/",
                BASE, STEP_INDEX))
