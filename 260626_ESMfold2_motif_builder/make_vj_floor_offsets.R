# =============================================================================
# Build the V/J FLOOR offset table read by ESM_motif_builder.R (VJ_FLOOR_TABLE).
#
# For each V gene, over a panel of BACKGROUND epitopes that share one step-0
# panel, take the gene's share of the pass set in every epitope and find the
# FLOOR, its smallest share across them. The floor is what the gene takes
# regardless of peptide. The offset is the score shift that would bring the
# floor epitope's share down to the uniform 1/n_genes baseline. Genes whose
# floor is already at or below baseline get no offset.
#
# Calibrating at the floor is the safety property: the penalty is sized at the
# gene's WORST epitope, so whatever it earns above its own floor survives.
#
# IMPORTANT - PASS must match how deep you actually cut. The offsets are only
# valid at the depth they were calibrated at: a gene's advantage grows as the
# cut deepens, so a table built at PASS = 0.20 over-penalises when the pipeline
# runs an absolute threshold that passes ~10%. Check the real depth with
#     quantile-of-threshold = mean(score >= ESM_THRESHOLDS[step])
# and set PASS to that.
#
# Usage:  Rscript make_vj_floor_offsets.R
# =============================================================================

# ---- config -----------------------------------------------------------------
BG_DIR    <- "VJ_background_correction"   # holds <epitope>/step0/{model,output}_<chain>.csv
BG_STEP   <- "step0"
EPITOPES  <- NULL                          # NULL = every A0201_* (or *_*) dir under BG_DIR
PASS      <- 0.10                          # pass fraction the offsets are calibrated at
SCORE_COL <- "iptm_pair_mean"
OUT_CSV   <- "vj_floor_offsets.csv"

# -----------------------------------------------------------------------------
if (is.null(EPITOPES)) {
  EPITOPES <- list.dirs(BG_DIR, recursive = FALSE, full.names = FALSE)
  EPITOPES <- EPITOPES[dir.exists(file.path(BG_DIR, EPITOPES, BG_STEP))]
}
if (length(EPITOPES) < 2)
  stop("Need at least 2 background epitopes under ", BG_DIR, "; found ", length(EPITOPES))

load_one <- function(e, chain) {
  d <- file.path(BG_DIR, e, BG_STEP)
  mf <- file.path(d, sprintf("model_%s.csv", chain))
  of <- file.path(d, sprintf("output_%s.csv", chain))
  if (!file.exists(mf) || !file.exists(of)) return(NULL)
  m <- read.csv(mf, stringsAsFactors = FALSE)
  o <- read.csv(of, stringsAsFactors = FALSE)
  names(o)[tolower(names(o)) == "id"] <- "id"
  x <- merge(m, o[, c("id", SCORE_COL)], by = "id")
  x$score <- x[[SCORE_COL]]
  # Both model tables carry BOTH chains, the partner suffixed "_dummy".
  x$V <- x[[if (chain == "alpha") "TRAV" else "TRBV"]]
  x <- x[!is.na(x$score) & !is.na(x$V), ]
  if (mean(grepl("_dummy$", x$V)) > 0.5)
    stop(sprintf("%s: the %s table's V column is the dummy chain.", e, chain))
  x
}

rows <- list()
for (chain in c("alpha", "beta")) {
  dat <- lapply(EPITOPES, load_one, chain = chain)
  names(dat) <- EPITOPES
  dat <- dat[!vapply(dat, is.null, logical(1))]
  if (length(dat) < 2) { warning("Skipping ", chain, ": fewer than 2 folded epitopes."); next }

  # panel identity: the floor is only meaningful on a shared panel
  keys <- lapply(dat, function(d) sort(d$id))
  if (!all(vapply(keys[-1], identical, logical(1), keys[[1]])))
    warning(sprintf("[%s] background epitopes do NOT share one step-0 panel - ",
                    chain), "floors are not strictly comparable.")

  genes <- sort(unique(dat[[1]]$V))
  base  <- 100 / length(genes)

  # share of the pass set, per gene per epitope
  P <- sapply(dat, function(d) {
    thr <- quantile(d$score, 1 - PASS, na.rm = TRUE)
    sel <- d$V[d$score >= thr]
    100 * as.numeric(table(factor(sel, levels = genes))) / length(sel)
  })
  rownames(P) <- genes

  floor_pct <- apply(P, 1, min)
  floor_ep  <- apply(P, 1, which.min)

  for (g in genes) {
    if (floor_pct[g] <= base) next
    d     <- dat[[floor_ep[g]]]
    thr   <- quantile(d$score, 1 - PASS, na.rm = TRUE)
    n_tar <- max(1, floor(base / 100 * sum(d$score >= thr)))   # folds it may keep
    v     <- sort(d$score[d$V == g], decreasing = TRUE)
    if (length(v) < n_tar) next
    delta <- max(0, v[n_tar] - thr)                            # shift to leave exactly n_tar
    if (delta <= 0) next
    rows[[length(rows) + 1]] <- data.frame(
      chain = chain, gene = g, delta = round(delta, 5),
      floor_pct = round(floor_pct[g], 2),
      floor_epitope = names(dat)[floor_ep[g]],
      baseline_pct = round(base, 2),
      max_pct = round(max(P[g, ]), 2),
      pass_calibrated = PASS, n_epitopes = length(dat),
      stringsAsFactors = FALSE)
  }
}

if (length(rows) == 0) stop("No gene had a floor above baseline - nothing to write.")
out <- do.call(rbind, rows)
write.csv(out, OUT_CSV, row.names = FALSE)

message(sprintf("PASS = %.2f, %d background epitopes: %s",
                PASS, length(EPITOPES), paste(EPITOPES, collapse = ", ")))
print(out, row.names = FALSE)
message(sprintf("\nWritten to: %s", OUT_CSV))
message("Set VJ_FLOOR_TABLE in ESM_motif_builder.R to this file.")
