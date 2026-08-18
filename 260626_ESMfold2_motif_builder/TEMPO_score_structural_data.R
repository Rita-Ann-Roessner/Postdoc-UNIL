# =============================================================================
# Score the STRUCTURAL motif(s) on the SAME per-fold test sets used by the
# experimental 5-fold CV — so structural vs experimental is a matched comparison
# (identical test rows, identical mean±band aggregation).
#
# Test folds are rebuilt from the AF3 benchmark table with the exact same
# deterministic filter the experimental script uses (CV_id == k & data_usage ==
# "test"), so the two are identical by construction — this script does NOT depend
# on the experimental script having run.
#
# The structural motif is FIXED (never trains on experimental TCRs), so per fold
# we just re-score it on that fold's test set. One structural motif per
# (epitope, step); path mirrors the notebook's atlas logic.
#
# SCORING ONLY — no AUC (the notebook draws ROC / computes AUC).
# Output: exp_data/<epitope>/struct_<step>_fold<k>_pred.csv   (Label + score)
# =============================================================================

suppressMessages(library(TEMPOtrain))

# ---- configuration ----------------------------------------------------------
#SUP_TABLE <- "../260330_af3_benchmark/Supplementary_Tables/Supplementary_Table_1.csv"
#EPITOPES  <- c("A0201_ELAGIGILTV", "A0201_GILGFVFTL", "A0201_LLWNGPMAV")

SUP_TABLE <- "IMMREP25/validation_all_cv.csv"
EPITOPES  <- c("A0201_KLYPFLWFA", "A0201_SQFNWTIYL", "A0201_TVYPYGTSL")
TOPDIR <- "IMMREP25"                 # holds <ep>/<step>/motif_top_binders.csv

OUT_DIR   <- "exp_data"
STEPS     <- c("step0", "step1", "step2")

# Chain combination — match the experimental script (and the structural pipeline).
CHAIN_WEIGHT <- "linear"

TCR_COLS <- c("TRAV", "TRAJ", "cdr3_TRA", "TRBV", "TRBJ", "cdr3_TRB")

# Structural motif path per (epitope, step).
motif_path <- function(ep, step) file.path(TOPDIR, ep, step, "motif_top_binders.csv")

# ---- TEMPO scoring: motif -> score test -> Label + score --------------------
score_fold <- function(motif_file, test_file, tmp_dir, chain_weight = CHAIN_WEIGHT) {
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  sc <- function(ch) {
    TEMPOtrain(input.train = motif_file, input.pred = test_file, output.path = tmp_dir,
               chain = ch, build.prank = FALSE, compute.prank = FALSE,
               write.data.pred = FALSE, write.data.train = FALSE, write.predictor = FALSE,
               return.object = TRUE)$score
  }
  if (chain_weight %in% c("sqrt", "linear")) {
    m  <- read.csv(motif_file)
    nA <- sum(!is.na(m$cdr3_TRA)); nB <- sum(!is.na(m$cdr3_TRB))
    wf <- if (chain_weight == "sqrt") sqrt else identity
    wA <- wf(nA); wB <- wf(nB)
    dA <- sc("A"); dB <- sc("B")
    if (nrow(dA) != nrow(dB)) stop("per-chain row counts differ (A=", nrow(dA), ", B=", nrow(dB), ")")
    pred <- dA; pred$score <- wA * dA$score + wB * dB$score
  } else {
    pred <- sc("AB")
  }
  pred
}

# ---- per epitope x fold x step ----------------------------------------------
d <- read.csv(SUP_TABLE, check.names = FALSE, stringsAsFactors = FALSE)

for (ep in EPITOPES) {
  g <- d[d$Epitope == ep, ]
  if (nrow(g) == 0) { warning("no rows for ", ep); next }
  ep_dir <- file.path(OUT_DIR, ep)
  dir.create(ep_dir, showWarnings = FALSE, recursive = TRUE)

  for (k in sort(unique(g$CV_id))) {
    te <- g[g$CV_id == k & g$data_usage == "test", ]                    # same filter as exp script
    if (nrow(te) == 0 || length(unique(te$Label)) < 2) {
      warning(sprintf("[%s] fold %s: skipped (test=%d, test_classes=%d)",
                      ep, k, nrow(te), length(unique(te$Label)))); next
    }
    fold_dir <- file.path(ep_dir, sprintf("fold%s", k))
    dir.create(fold_dir, showWarnings = FALSE, recursive = TRUE)

    test <- te[, c(TCR_COLS, "Label")]; test$model <- ep               # identical to exp fold test
    tf <- file.path(fold_dir, "test.csv"); write.csv(test, tf, row.names = FALSE)

    for (step in STEPS) {
      mfile <- motif_path(ep, step)
      if (!file.exists(mfile)) { warning(sprintf("[%s|%s] motif missing: %s", ep, step, mfile)); next }
      pred <- score_fold(mfile, tf, file.path(fold_dir, sprintf("struct_%s_tmp", step)))
      write.csv(pred[, c("Label", "score")],
                file.path(ep_dir, sprintf("struct_%s_fold%s_pred.csv", step, k)), row.names = FALSE)
    }
    message(sprintf("[%s] fold %s: test=%d (pos=%d) scored for %d step(s)",
                    ep, k, nrow(te), sum(te$Label == 1), length(STEPS)))
  }
}

message("\nDone. Per-fold structural predictions in ",
        OUT_DIR, "/<epitope>/struct_<step>_fold<k>_pred.csv")
