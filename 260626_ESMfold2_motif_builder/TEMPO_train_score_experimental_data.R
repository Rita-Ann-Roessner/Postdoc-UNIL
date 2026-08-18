# =============================================================================
# Train + score TEMPO on the EXPERIMENTAL reference TCRs with 5-fold CV.
#
# Leakage-free split comes straight from the AF3 benchmark table:
#   data_usage  -> "train" / "test"   (which rows train the motif vs are scored)
#   CV_id       -> 1..5                (the cross-validation fold)
# For each fold k we restrict to CV_id == k and, WITHIN that fold, train on the
# `train` positives and score the `test` rows. Within a CV_id the train and test
# TCRs are clonotype-disjoint (verified), so there is no train/test leakage.
#
# SCORING ONLY — no AUC (the notebook draws ROC / computes AUC).
# Output: exp_data/<epitope>/exp_fold<k>_pred.csv   (Label + score, one per fold)
# (per-fold train/test/temp files go in exp_data/<epitope>/fold<k>/)
# =============================================================================

suppressMessages(library(TEMPOtrain))

# ---- configuration ----------------------------------------------------------
#SUP_TABLE <- "../260330_af3_benchmark/Supplementary_Tables/Supplementary_Table_1.csv"
#EPITOPES  <- c("A0201_ELAGIGILTV", "A0201_GILGFVFTL", "A0201_LLWNGPMAV")

SUP_TABLE <- "IMMREP25/validation_all_cv.csv"
EPITOPES  <- c("A0201_KLYPFLWFA", "A0201_SQFNWTIYL", "A0201_TVYPYGTSL")

OUT_DIR   <- "exp_data"


# Chain combination — "linear"/"sqrt" weight each chain by motif passer count
# (as in the structural pipeline); "equal" = plain TEMPO AB score.
CHAIN_WEIGHT <- "linear"

TCR_COLS <- c("TRAV", "TRAJ", "cdr3_TRA", "TRBV", "TRBJ", "cdr3_TRB")

# ---- TEMPO scoring for one fold: train motif -> score test -> Label + score ---
score_fold <- function(train_motif, test_file, tmp_dir, chain_weight = CHAIN_WEIGHT) {
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  sc <- function(ch) {
    TEMPOtrain(input.train = train_motif, input.pred = test_file, output.path = tmp_dir,
               chain = ch, build.prank = FALSE, compute.prank = FALSE,
               write.data.pred = FALSE, write.data.train = FALSE, write.predictor = FALSE,
               return.object = TRUE)$score
  }
  if (chain_weight %in% c("sqrt", "linear")) {
    m  <- read.csv(train_motif)
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

# ---- per epitope x fold -----------------------------------------------------
d <- read.csv(SUP_TABLE, check.names = FALSE, stringsAsFactors = FALSE)

for (ep in EPITOPES) {
  g <- d[d$Epitope == ep, ]
  if (nrow(g) == 0) { warning("no rows for ", ep); next }
  ep_dir <- file.path(OUT_DIR, ep)
  dir.create(ep_dir, showWarnings = FALSE, recursive = TRUE)

  for (k in sort(unique(g$CV_id))) {
    tr <- g[g$CV_id == k & g$data_usage == "train" & g$Label == 1, ]   # motif = train positives
    te <- g[g$CV_id == k & g$data_usage == "test", ]                    # test = both classes
    if (nrow(tr) < 2 || nrow(te) == 0 || length(unique(te$Label)) < 2) {
      warning(sprintf("[%s] fold %s: skipped (train_pos=%d, test=%d, test_classes=%d)",
                      ep, k, nrow(tr), nrow(te), length(unique(te$Label)))); next
    }

    fold_dir <- file.path(ep_dir, sprintf("fold%s", k))
    dir.create(fold_dir, showWarnings = FALSE, recursive = TRUE)

    motif <- tr[, TCR_COLS]; motif$model <- ep                          # TEMPO input.train
    test  <- te[, c(TCR_COLS, "Label")]; test$model <- ep               # TEMPO input.pred
    mf <- file.path(fold_dir, "train_motif.csv"); write.csv(motif, mf, row.names = FALSE)
    tf <- file.path(fold_dir, "test.csv");        write.csv(test,  tf, row.names = FALSE)

    pred <- score_fold(mf, tf, fold_dir)
    write.csv(pred[, c("Label", "score")],
              file.path(ep_dir, sprintf("exp_fold%s_pred.csv", k)), row.names = FALSE)
    message(sprintf("[%s] fold %s: train_pos=%d  test=%d (pos=%d)",
                    ep, k, nrow(tr), nrow(te), sum(te$Label == 1)))
  }
}

message("\nDone. Per-fold predictions in ", OUT_DIR, "/<epitope>/exp_fold<k>_pred.csv")
