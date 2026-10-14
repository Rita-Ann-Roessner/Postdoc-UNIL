library(TEMPOtrain)
library(MixTCRviz)
library(pROC)

# =============================================================================
# Motif builder benchmark — step-0 analysis
#
# Analysis only: no TCR generation, no input building for later steps. For every
# method and epitope it reads the existing step-0 fold scores and writes
#   step0/top_tcrs_<chain>.csv   — the passing TCRs per chain
#   step0/motif_<chain>/         — MixTCRviz motif of those passers
#   step0/motif_top_binders.csv  — combined TEMPO training set
#   step0/validation/validation_pred.csv  — scored validation set (+ AUC)
# plus one benchmark_auc_summary.csv at the top level.
#
# Expected layout (already produced by the generation pipeline):
#   <method>/<MHC_PEPTIDE>/validation.csv
#   <method>/<MHC_PEPTIDE>/step0/{input,model}_<chain>.csv, output_<chain>.csv
#
# Usage: set RUN_METHODS in the configuration block below, then source the script
# (RStudio: Source, or Ctrl/Cmd-Shift-S). The working directory must be the one
# holding the <method>/ directories — in RStudio, Session > Set Working Directory
# > To Source File Location.
# =============================================================================

### ---- Configuration --------------------------------------------------------

STEP_DIR  <- "step0"            # sub-directory holding the folds to analyse
SCORE_COL <- "iptm_pair_mean"   # column in output_<chain>.csv; higher = better

# Selection is by percentile, not by absolute score: the score scales are not
# comparable across folding engines (step-0 median iptm_pair_mean runs ~0.11 for
# AF3 vs ~0.57 for OpenDDE), so a shared absolute cutoff would select 5% of one
# method's TCRs and 62% of another's, and the AUCs would differ by pass rate
# rather than by motif quality. A shared percentile selects the same FRACTION of
# each method's own distribution, so every method contributes an equally sized
# motif and only the ranking of its TCRs is being compared.
#
# PCTILE = 84 reproduces the ESMFold2 reference cutoff of 0.50 on the alpha
# chain (cut = 0.501, 363 vs 368 passers). Its beta chain matches 0.50 at p76.6
# instead (beta scores sit lower relative to 0.5), so the two chains cannot both
# be matched by one value; alpha is the one tracked here, since it dominates both
# the motif (368 of 505 rows) and the "linear" chain weighting. The pooled
# alpha+beta match is p82.3.
PCTILE <- 84

# One entry per method, each taking EITHER
#   percentile = keep the top (100 - p)% of that chain's own score distribution
#                (computed per chain, on that chain's own scores), OR
#   threshold  = fixed score cutoff (TCRs with score >= cutoff pass).
# percentile takes precedence when both are set.
METHODS <- list(
  ESMFold2 = list(dir = "ESMFold2", percentile = PCTILE),
  AF3      = list(dir = "AF3",      percentile = PCTILE),
  OpenDDE  = list(dir = "OpenDDE",  percentile = PCTILE),
  Protenix = list(dir = "Protenix", percentile = PCTILE)
  # Fixed-cutoff alternative, e.g.: list(dir = "AF3", threshold = 0.30)
)

# Which of the above to actually run — edit this line to run one method at a
# time. Any subset of names(METHODS); "all" (or NULL) runs every one.
#   RUN_METHODS <- "ESMFold2"               # one method
#   RUN_METHODS <- c("AF3", "Protenix")     # a subset
#   RUN_METHODS <- "all"                    # the full benchmark
RUN_METHODS <- "ESMFold2"

EPITOPES   <- NULL   # NULL = every sub-directory of each method dir
CHAINS     <- c(alpha = "A", beta = "B")

MIN_TOP_TCRS <- 10    # below this a chain is skipped (no motif, no training rows)
PLOT_MOTIF   <- TRUE  # write the MixTCRviz motif_<chain>/ plots
VALIDATE     <- TRUE  # run TEMPO validation (needs <method>/<epitope>/validation.csv)

CHAIN_WEIGHT <- "linear"  # how TEMPO combines the chains when scoring validation:
                          # "equal"  = default TEMPO AB score (equal weight).
                          # "sqrt"   = weight each chain by sqrt(n_passers).
                          # "linear" = weight by the passer count itself.

SUMMARY_FILE <- "benchmark_auc_summary.csv"


# =============================================================================
# Score loading
# =============================================================================

load_scores <- function(scores_file, model_file, score_col = SCORE_COL) {
  scores <- read.csv(scores_file, stringsAsFactors = FALSE, check.names = FALSE)
  model  <- read.csv(model_file,  stringsAsFactors = FALSE, check.names = FALSE)

  # Fold output keys on 'ID', the model table on 'id' — normalise both.
  fix_id <- function(df) {
    hit <- which(tolower(names(df)) == "id")
    if (length(hit) >= 1) names(df)[hit[1]] <- "id"
    df
  }
  scores <- fix_id(scores)
  model  <- fix_id(model)

  if (!score_col %in% colnames(scores))
    stop("Column '", score_col, "' not found in ", scores_file,
         "\nAvailable: ", paste(colnames(scores), collapse = ", "))
  if (!"id" %in% colnames(scores) || !"id" %in% colnames(model))
    stop("Both the scores file and the model table need an 'id'/'ID' column.")

  merge(model, scores[, c("id", score_col)], by = "id")
}


GENE_COLS <- c("TRAV", "TRAJ", "cdr3_TRA", "TRBV", "TRBJ", "cdr3_TRB")

# Metadata table matching a folded output_<chain>.csv. Prefer the immutable
# input_<chain>.csv (the exact batch that was folded) so scores are never joined
# against a model_<chain>.csv a later run may have regenerated — but only when it
# really is the gene table: for some methods input_<chain>.csv is the fold.py
# sequence file (ID, TCRA, TCRB, MHC, B2M, PEPTIDE), and then model_<chain>.csv
# is the only table carrying V/J + CDR3 (its ids match the folded batch).
meta_csv <- function(step_dir, chain_name) {
  inp <- file.path(step_dir, sprintf("input_%s.csv", chain_name))
  if (file.exists(inp)) {
    hdr <- names(read.csv(inp, nrows = 1, check.names = FALSE))
    if (all(GENE_COLS %in% hdr)) return(inp)
  }
  file.path(step_dir, sprintf("model_%s.csv", chain_name))
}


# Score cutoff for one chain, from this method's threshold / percentile setting.
chain_cutoff <- function(scored, cfg, label) {
  if (!is.null(cfg$percentile)) {
    cut <- as.numeric(quantile(scored[[SCORE_COL]], cfg$percentile / 100, na.rm = TRUE))
    message(sprintf("    cutoff = %.3f (%gth pctile of %d scores)",
                    cut, cfg$percentile, nrow(scored)))
    return(cut)
  }
  if (is.null(cfg$threshold))
    stop("Method config for '", label, "' sets neither 'threshold' nor 'percentile'.")
  cfg$threshold
}


# =============================================================================
# Motif helpers
# =============================================================================

# MixTCRviz motif for one chain's passers vs the default baseline.
plot_motif <- function(top_tcrs, chain_name, label, step_dir) {
  out_dir <- file.path(step_dir, sprintf("motif_%s", chain_name))
  message(sprintf("    plotting motif_%s (n=%d)", chain_name, nrow(top_tcrs)))
  MixTCRviz(
    input1      = top_tcrs,
    output.path = out_dir,
    plot        = TRUE,
    renormVJ    = TRUE,
    set.title   = sprintf("%s - %s top binders (n=%d)", label, chain_name, nrow(top_tcrs)),
    verbose     = 0
  )
}


# Combine alpha and beta passers into one TEMPO training df, clearing the
# dummy-chain columns in each (each row carries only its real chain).
build_motif_df <- function(top_alpha, top_beta, model_label) {
  cols <- GENE_COLS
  rows <- list()
  if (!is.null(top_alpha) && nrow(top_alpha) > 0) {
    df_a <- top_alpha[, cols, drop = FALSE]
    df_a$TRBV <- df_a$TRBJ <- df_a$cdr3_TRB <- NA_character_
    df_a$model <- model_label
    rows[["alpha"]] <- df_a
  }
  if (!is.null(top_beta) && nrow(top_beta) > 0) {
    df_b <- top_beta[, cols, drop = FALSE]
    df_b$TRAV <- df_b$TRAJ <- df_b$cdr3_TRA <- NA_character_
    df_b$model <- model_label
    rows[["beta"]] <- df_b
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}


# =============================================================================
# TEMPO validation
# =============================================================================

validate_motif <- function(motif_file, validation_file, output_dir,
                            peptide, mhc, score_col = "score", plot = TRUE,
                            chain_weight = CHAIN_WEIGHT) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Use TEMPO's returned data frame directly (res$score with Label + score) —
  # no disk glob, so no stale-file / locale-sort ambiguity.
  score_chain <- function(ch) {
    res <- TEMPOtrain(
      input.train = motif_file, input.pred = validation_file, output.path = output_dir,
      chain = ch, build.prank = FALSE, compute.prank = FALSE,
      write.data.pred = FALSE, write.data.train = FALSE, write.predictor = FALSE,
      return.object = TRUE)
    res$score
  }

  if (chain_weight %in% c("sqrt", "linear")) {
    # score = wA*score_A + wB*score_B, with w = n (linear) or sqrt(n) (gentler).
    motif <- read.csv(motif_file)
    nA <- sum(!is.na(motif$cdr3_TRA)); nB <- sum(!is.na(motif$cdr3_TRB))
    wfun <- if (chain_weight == "sqrt") sqrt else identity
    wA <- wfun(nA); wB <- wfun(nB)
    dA <- score_chain("A"); dB <- score_chain("B")
    if (nrow(dA) != nrow(dB))
      stop("Per-chain prediction row counts differ (A=", nrow(dA), ", B=", nrow(dB), ").")
    pred_df <- dA
    pred_df[[score_col]] <- wA * dA[[score_col]] + wB * dB[[score_col]]
    message(sprintf("    chain weights (%s): A=%.2f (n=%d)  B=%.2f (n=%d)",
                    chain_weight, wA, nA, wB, nB))
  } else {
    pred_df <- score_chain("AB")   # equal weight = default combined AB score
  }

  write.csv(pred_df, file.path(output_dir, "validation_pred.csv"), row.names = FALSE)

  if (!"Label" %in% colnames(pred_df))
    stop("'Label' column not found in validation output; check validation.csv.")
  if (!score_col %in% colnames(pred_df))
    stop("Column '", score_col, "' not found in validation output.\n",
         "Available: ", paste(colnames(pred_df), collapse = ", "))

  pred_df   <- pred_df[!is.na(pred_df[[score_col]]), ]
  roc_obj   <- roc(pred_df$Label, pred_df[[score_col]], quiet = TRUE)
  auc_val   <- as.numeric(auc(roc_obj))
  auc01_val <- as.numeric(auc(roc_obj,
                               partial.auc         = c(1, 0.9),
                               partial.auc.correct = TRUE,
                               partial.auc.focus   = "specificity"))

  message(sprintf("    AUC = %.4f  |  AUC0.1 = %.4f", auc_val, auc01_val))

  if (plot) {
    pdf(file.path(output_dir, "roc.pdf"), width = 5, height = 5)
    plot(roc_obj,
         main = sprintf("ROC - %s_%s", mhc, peptide),
         col  = "#2166ac", lwd = 2, xaxt = "n",
         xlab = "FPR (1 - Specificity)", ylab = "TPR (Sensitivity)")
    axis(1, at = seq(1, 0, -0.2), labels = seq(0, 1, 0.2))
    legend("bottomright", bty = "o", bg = "white", box.col = "grey80",
           legend = c(sprintf("AUC    = %.4f", auc_val),
                      sprintf("AUC0.1 = %.4f", auc01_val)))
    dev.off()
  }

  list(auc = auc_val, auc01 = auc01_val)
}


# =============================================================================
# Per-epitope analysis
# =============================================================================

analyse_epitope <- function(method, cfg, epitope) {
  mhc      <- sub("_.*", "", epitope)
  peptide  <- sub("^[^_]*_", "", epitope)
  step_dir <- file.path(cfg$dir, epitope, STEP_DIR)

  message(sprintf("  [%s / %s]", method, epitope))
  if (!dir.exists(step_dir)) {
    warning(sprintf("[%s/%s] No %s directory — skipped.", method, epitope, STEP_DIR))
    return(NULL)
  }

  # ---- per-chain filtering ---------------------------------------------------
  tops <- list()
  for (chain_name in names(CHAINS)) {
    scores_file <- file.path(step_dir, sprintf("output_%s.csv", chain_name))
    if (!file.exists(scores_file)) {
      warning(sprintf("[%s/%s] Missing %s — chain skipped.", method, epitope, scores_file))
      next
    }
    scored   <- load_scores(scores_file, meta_csv(step_dir, chain_name))
    cutoff   <- chain_cutoff(scored, cfg, method)
    top_tcrs <- scored[scored[[SCORE_COL]] >= cutoff, ]
    message(sprintf("    %-5s %d / %d TCRs with %s >= %.3f",
                    chain_name, nrow(top_tcrs), nrow(scored), SCORE_COL, cutoff))

    write.csv(top_tcrs, file.path(step_dir, sprintf("top_tcrs_%s.csv", chain_name)),
              row.names = FALSE)

    if (nrow(top_tcrs) < MIN_TOP_TCRS) {
      warning(sprintf("[%s/%s] Only %d passers for %s (< %d) — chain dropped.",
                      method, epitope, nrow(top_tcrs), chain_name, MIN_TOP_TCRS))
      next
    }
    tops[[chain_name]] <- top_tcrs
    if (PLOT_MOTIF) plot_motif(top_tcrs, chain_name, epitope, step_dir)
  }

  # ---- combined motif training set -------------------------------------------
  motif_df <- build_motif_df(tops$alpha, tops$beta, epitope)
  if (is.null(motif_df)) {
    warning(sprintf("[%s/%s] No passing TCRs on either chain — nothing written.", method, epitope))
    return(NULL)
  }
  motif_file <- file.path(step_dir, "motif_top_binders.csv")
  write.csv(motif_df, motif_file, row.names = FALSE)
  message(sprintf("    motif_top_binders.csv: %d TCRs", nrow(motif_df)))

  # ---- validation ------------------------------------------------------------
  if (!VALIDATE) return(NULL)
  validation_file <- file.path(cfg$dir, epitope, "validation.csv")
  if (!file.exists(validation_file)) {
    warning(sprintf("[%s/%s] No validation.csv — motif written, validation skipped.",
                    method, epitope))
    return(NULL)
  }
  res <- validate_motif(
    motif_file      = motif_file,
    validation_file = validation_file,
    output_dir      = file.path(step_dir, "validation"),
    peptide         = peptide,
    mhc             = mhc
  )
  data.frame(method = method, epitope = epitope,
             n_motif = nrow(motif_df),
             n_alpha = if (is.null(tops$alpha)) 0L else nrow(tops$alpha),
             n_beta  = if (is.null(tops$beta))  0L else nrow(tops$beta),
             auc = res$auc, auc01 = res$auc01, row.names = NULL)
}


# =============================================================================
# Main
# =============================================================================

selected <- RUN_METHODS
if (is.null(selected) || identical(tolower(selected), "all")) {
  selected <- names(METHODS)
}

unknown <- setdiff(selected, names(METHODS))
if (length(unknown) > 0)
  stop("Unknown method(s): ", paste(unknown, collapse = ", "),
       "\nAvailable: ", paste(names(METHODS), collapse = ", "))

# A partial run must not overwrite the full-benchmark summary.
summary_file <- if (length(selected) == length(METHODS)) {
  SUMMARY_FILE
} else {
  sub("\\.csv$", sprintf("_%s.csv", paste(selected, collapse = "-")), SUMMARY_FILE)
}

message(sprintf("Running method(s): %s  ->  %s",
                paste(selected, collapse = ", "), summary_file))

summary_rows <- list()

for (method in selected) {
  cfg <- METHODS[[method]]
  if (!dir.exists(cfg$dir)) {
    warning(sprintf("Method directory missing: %s — skipped.", cfg$dir))
    next
  }
  eps <- if (is.null(EPITOPES)) {
    sort(basename(list.dirs(cfg$dir, recursive = FALSE)))
  } else EPITOPES

  message(sprintf("\n========== %s (%d epitope%s) ==========",
                  method, length(eps), if (length(eps) == 1) "" else "s"))

  for (epitope in eps) {
    row <- analyse_epitope(method, cfg, epitope)
    if (!is.null(row)) summary_rows[[paste(method, epitope)]] <- row
  }
}

if (length(summary_rows) > 0) {
  auc_df <- do.call(rbind, summary_rows)
  write.csv(auc_df, summary_file, row.names = FALSE)
  message("\n===== AUC summary =====")
  for (i in seq_len(nrow(auc_df)))
    message(sprintf("  %-10s %-18s AUC = %.4f  |  AUC0.1 = %.4f  (n_motif = %d)",
                    auc_df$method[i], auc_df$epitope[i],
                    auc_df$auc[i], auc_df$auc01[i], auc_df$n_motif[i]))
  message(sprintf("  Saved to: %s", summary_file))
} else {
  message("\nNo validated epitopes — nothing written to the summary.")
}
