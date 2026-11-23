# =============================================================================
# score_tcrs.R — score an arbitrary TCR CSV against a trained (STEP="final") motif
# =============================================================================
# Prereq: run ESM_motif_builder.R with STEP <- "final" for the epitope first, so
#         the motif file <BASE_OUTPUT_DIR>/<label>/step<N_STEPS>/motif_top_binders.csv
#         exists.
#
# Input CSV must contain (any order): TRAV, cdr3_TRA, TRBV, cdr3_TRB, TRAJ, TRBJ.
# A 'model' column is not required; one is added if absent. An existing 'model'
# column is passed through untouched and need NOT match the motif's own model
# label — TEMPO does not use it to match prediction rows to the trained model
# (verified: a query labelled "LAU155" scored normally against the
# "A0201_SLLMWITQC" motif). It is dropped from the output, see DROP_MODEL_COL.
# A 'Label' column is not required either; nothing here computes an AUC.
#
# Edit the CONFIG block below, then source / run the whole script in RStudio.
# =============================================================================

# ============================ CONFIG =========================================
TCR_CSV      <- "/Users/roessner/Documents/PostDoc/Data/260626_ESMfold2_motif_builder/TCR_motif_atlas_cdr3_templated_new_installation_reduced_cycles_NY-ESO-1_fix/LAU50.csv"                   # TCRs to score (6 gene/CDR3 cols)
MOTIF_FILE   <- "/Users/roessner/Documents/PostDoc/Data/260626_ESMfold2_motif_builder/TCR_motif_atlas_cdr3_templated_new_installation_reduced_cycles_NY-ESO-1_fix/A0201_SLLMWITQC/step2/motif_top_binders.csv"
OUT_CSV      <- "/Users/roessner/Documents/PostDoc/Data/260626_ESMfold2_motif_builder/TCR_motif_atlas_cdr3_templated_new_installation_reduced_cycles_NY-ESO-1_fix/LAU50_A0201_SLLMWITQC_model.csv"                            # NULL = "<TCR_CSV>_scored.csv"
CHAIN_WEIGHT <- "linear"                        # equal | sqrt | linear — MUST match the
                                                # CHAIN_WEIGHT the motif was validated with
                                                # in ESM_motif_builder.R (not read from it)
DROP_MODEL_COL <- TRUE                          # drop the passed-through 'model' column from
                                                # the output (it is provenance of the INPUT,
                                                # not of the model that did the scoring)

BUILD_PRANK  <- TRUE     # ask TEMPO for percentile ranks (build.prank + compute.prank).
                         # build.prank builds TEMPO's own rank reference (understood to be its
                         # ~1M-TCR background repertoire, which is why it is separate from
                         # compute.prank: expensive to build, cheap to apply, and savable with
                         # the predictor). The ranks are PER CHAIN, on the raw per-chain score.
                         # Sanity check on the first run: 21 query TCRs ranked against a large
                         # background should give clustered, arbitrary percentages. If instead
                         # they come back as near-even multiples of 100/21 (0, 4.8, 9.5, ...),
                         # the reference is the query itself and not a repertoire.

BACKGROUND_CSV <- NULL   # optional repertoire CSV (same 6 gene/CDR3 columns). This does NOT
                         # replace TEMPO's prank; it adds 'prank_weighted' for the CHAIN_WEIGHT
                         # combined score, which TEMPO structurally cannot rank: the weighted
                         # score is formed here, after both per-chain calls have returned, so
                         # TEMPO never sees it. Leave NULL if the per-chain pranks are enough.
# =============================================================================

# Pull the scored data frame out of TEMPO's returned object, plus any percentile-rank
# columns it carries. The prank column name is discovered rather than assumed, and
# reported once, so a naming change in TEMPOtrain surfaces instead of going missing.
.tempo_frame <- function(res, chain, tmp_dir = NULL, written = character(0)) {
  d <- res$score
  if (!is.data.frame(d)) stop("TEMPOtrain's res$score is not a data frame (chain ", chain, ").")
  # TEMPO's own name for it is 'perc_rank'. Kept broad so a rename still matches.
  pfind <- function(x) grep("perc.?rank|prank|percentile", x, ignore.case = TRUE, value = TRUE)

  pr <- pfind(colnames(d))
  src <- "res$score"

  # 2nd: a same-length numeric sibling on the returned object.
  if (!length(pr)) {
    cand <- pfind(names(res)[vapply(res, function(x)
      is.numeric(x) && length(x) == nrow(d), logical(1))])
    for (nm in cand) d[[nm]] <- res[[nm]]
    pr <- cand; src <- "res$<sibling>"
  }

  # 3rd: the prediction table TEMPO wrote to disk. Joined on .row, never on
  # position, so a reordered or partial file cannot silently misalign.
  if (!length(pr) && !is.null(tmp_dir)) {
    new_files <- setdiff(list.files(tmp_dir, pattern = "\\.(csv|txt|tsv)$",
                                    full.names = TRUE, recursive = TRUE), written)
    for (f in new_files) {
      w <- tryCatch(read.csv(f, stringsAsFactors = FALSE, check.names = FALSE),
                    error = function(e) NULL)
      if (is.null(w) || !".row" %in% colnames(w)) next
      got <- pfind(colnames(w))
      if (!length(got)) next
      idx <- match(d$.row, w$.row)
      if (anyNA(idx)) next
      for (nm in got) d[[nm]] <- w[[nm]][idx]
      pr <- got; src <- basename(f)
      break
    }
  }

  message(sprintf("  [chain %s] prank: %s", chain,
                  if (length(pr)) sprintf("%s (from %s)", paste(pr, collapse = ", "), src)
                  else "NONE found in the returned object or the written files"))
  attr(d, "prank_cols") <- pr
  d
}

score_tcrs <- function(tcr_csv, motif_file,
                       out_csv        = NULL,
                       chain_weight   = "linear",
                       drop_model_col = TRUE,
                       build_prank    = TRUE,
                       background_csv = NULL,
                       tmp_dir        = tempfile("score_tcrs_")) {

  suppressMessages(library(TEMPOtrain))
  stopifnot(file.exists(tcr_csv), file.exists(motif_file))
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

  req <- c("TRAV", "TRAJ", "cdr3_TRA", "TRBV", "TRBJ", "cdr3_TRB")
  tcrs <- read.csv(tcr_csv, stringsAsFactors = FALSE, check.names = FALSE)
  miss <- setdiff(req, colnames(tcrs))
  if (length(miss)) stop("Input CSV missing columns: ", paste(miss, collapse = ", "))

  # TEMPO's input.pred needs a 'model' column; keep original row order via an index.
  # .row is load-bearing: TEMPO's returned frame is NOT assumed to preserve input
  # order, and .row is what restores it. It survives the round trip as an extra
  # passthrough column (verified), and is checked for below rather than trusted.
  pred <- tcrs
  pred$.row  <- seq_len(nrow(pred))
  if (!"model" %in% colnames(pred)) pred$model <- "query"
  pred_file <- file.path(tmp_dir, "pred.csv")
  write.csv(pred, pred_file, row.names = FALSE)

  # Weights come from the MOTIF (the training set), never from the query.
  motif <- read.csv(motif_file, stringsAsFactors = FALSE, check.names = FALSE)
  nA <- sum(!is.na(motif$cdr3_TRA)); nB <- sum(!is.na(motif$cdr3_TRB))
  wfun <- if (chain_weight == "sqrt") sqrt else identity
  wA <- if (chain_weight %in% c("sqrt", "linear")) wfun(nA) else 1
  wB <- if (chain_weight %in% c("sqrt", "linear")) wfun(nB) else 1

  score_chain <- function(ch, pred_path, prank) {
    # write.data.pred follows the prank request: if TEMPO only reports the rank in
    # the written prediction table (as the celiac script's all-TRUE call suggests
    # it might), that file is the fallback source. Written into tmp_dir, not the
    # project. Snapshot first so the fallback reads only this call's output.
    before <- list.files(tmp_dir, pattern = "\\.(csv|txt|tsv)$",
                         full.names = TRUE, recursive = TRUE)
    res <- TEMPOtrain(
      input.train = motif_file, input.pred = pred_path, output.path = tmp_dir,
      chain = ch, build.prank = prank, compute.prank = prank,
      filename.train = "motif", filename.pred = paste0("query_", ch),
      write.data.pred = isTRUE(prank), write.data.train = FALSE,
      write.predictor = FALSE, return.object = TRUE)
    .tempo_frame(res, ch, tmp_dir = if (isTRUE(prank)) tmp_dir else NULL,
                 written = before)
  }

  # Score one prediction CSV with the configured chain rule. Returns a frame with
  # 'score' (weighted where applicable) plus any per-chain prank columns.
  score_set <- function(pred_path, prank) {
    if (chain_weight %in% c("sqrt", "linear")) {
      dA <- score_chain("A", pred_path, prank)
      dB <- score_chain("B", pred_path, prank)
      if (nrow(dA) != nrow(dB))
        stop("Per-chain row counts differ (A=", nrow(dA), ", B=", nrow(dB), ").")
      # Suffix EVERY per-chain column, not just the ones a pattern recognises.
      # Anything TEMPO returns (score, perc_rank, Problem, ...) is a property of the
      # chain it was computed for, so an unsuffixed copy would be chain A's value
      # masquerading as the TCR's, with chain B's silently dropped.
      keys <- intersect(c(".row", "Label", "model"), colnames(dA))
      o <- dA[, keys, drop = FALSE]
      for (nm in setdiff(colnames(dA), keys)) o[[paste0(nm, "_A")]] <- dA[[nm]]
      for (nm in setdiff(colnames(dB), keys)) o[[paste0(nm, "_B")]] <- dB[[nm]]
      o$score <- wA * dA$score + wB * dB$score   # the combined score, after the parts
      o
    } else {
      score_chain("AB", pred_path, prank)   # equal-weight combined score
    }
  }

  if (chain_weight %in% c("sqrt", "linear")) {
    message(sprintf("chain weights (%s): A=%.2f (n=%d)  B=%.2f (n=%d)",
                    chain_weight, wA, nA, wB, nB))
    if (isTRUE(all.equal(wA, wB)))
      message("  note: wA == wB — the weighting does not change the RANKING here.")
  }

  message("Scoring query (", nrow(tcrs), " TCRs):")
  out <- score_set(pred_file, build_prank)

  # Guard before aligning scores back onto the input. A row TEMPO silently drops
  # (unparseable CDR3, unknown gene) would otherwise shift every score by one and
  # still produce a plausible-looking table — these scores are large, unbounded
  # and have no reference value to eyeball against, so it would not be noticed.
  if (!".row" %in% colnames(out))
    stop("TEMPO's output lost the .row index — cannot verify row alignment.")
  if (nrow(out) != nrow(tcrs))
    stop("TEMPO returned ", nrow(out), " rows for ", nrow(tcrs), " input TCRs; ",
         "scores cannot be aligned. Check the input for unparseable genes/CDR3s.")

  # Restore original input order and attach the score to the original columns.
  out <- out[order(out$.row), , drop = FALSE]
  res <- tcrs
  res$score <- out$score

  # Per-chain scores and every prank column TEMPO returned, carried through.
  extra <- setdiff(colnames(out), c(colnames(tcrs), ".row", "score", "model", "Label"))
  for (nm in extra) res[[nm]] <- out[[nm]]

  # Is TEMPO's prank really against its own background, or against this query?
  # Ranks over a large reference are clustered and arbitrary; ranks computed over
  # the n query TCRs alone are n evenly spaced multiples of 100/n. Cheap to tell apart.
  for (nm in grep("perc.?rank|prank|percentile", extra, ignore.case = TRUE, value = TRUE)) {
    v <- sort(unique(res[[nm]][!is.na(res[[nm]])]))
    if (length(v) == nrow(res) && length(v) > 2) {
      gaps <- diff(v)
      if (stats::sd(gaps) < 0.01 * mean(gaps))
        warning("'", nm, "' is ", length(v), " evenly spaced values (step ",
                sprintf("%.2f", mean(gaps)), " ~ 100/", nrow(res),
                "): the rank reference looks like the QUERY, not a repertoire.")
    }
  }

  # Percentile rank of the WEIGHTED score against a background repertoire. Computed
  # here rather than by TEMPO, because the weighted score only exists after both
  # chains have come back, so TEMPO cannot rank it.
  if (!is.null(background_csv)) {
    stopifnot(file.exists(background_csv))
    bg <- read.csv(background_csv, stringsAsFactors = FALSE, check.names = FALSE)
    bmiss <- setdiff(req, colnames(bg))
    if (length(bmiss)) stop("BACKGROUND_CSV missing columns: ", paste(bmiss, collapse = ", "))
    bg$.row <- seq_len(nrow(bg))
    if (!"model" %in% colnames(bg)) bg$model <- "background"
    bg_file <- file.path(tmp_dir, "background.csv")
    write.csv(bg, bg_file, row.names = FALSE)

    message("Scoring background repertoire (", nrow(bg), " TCRs):")
    bg_out <- score_set(bg_file, prank = FALSE)   # only the score is needed as reference
    ref <- sort(bg_out$score[!is.na(bg_out$score)])
    if (length(ref) < 100)
      warning("background has only ", length(ref),
              " scorable TCRs — prank_weighted will be coarse.")
    # % of background scoring strictly below each query TCR.
    res$prank_weighted <- 100 * findInterval(res$score, ref) / length(ref)
    message(sprintf("  prank_weighted vs %d background TCRs (median query rank %.1f%%)",
                    length(ref), stats::median(res$prank_weighted, na.rm = TRUE)))
  }

  if (any(is.na(res$score)))
    warning(sum(is.na(res$score)), " of ", nrow(res),
            " TCRs scored NA (they sort last).")
  if (isTRUE(drop_model_col)) res$model <- NULL
  res <- res[order(-res$score), ]          # best binders first

  if (!is.null(out_csv)) {
    write.csv(res, out_csv, row.names = FALSE)
    message("Scored ", nrow(res), " TCRs -> ", out_csv)
  }
  res
}

# -----------------------------------------------------------------------------
# Run (uses the CONFIG block at the top)
# -----------------------------------------------------------------------------
out_csv <- if (is.null(OUT_CSV)) sub("\\.csv$", "_scored.csv", TCR_CSV) else OUT_CSV

scored <- score_tcrs(tcr_csv = TCR_CSV, motif_file = MOTIF_FILE,
                     out_csv = out_csv, chain_weight = CHAIN_WEIGHT,
                     drop_model_col = DROP_MODEL_COL,
                     build_prank = BUILD_PRANK, background_csv = BACKGROUND_CSV)
print(utils::head(scored, 20))
