# =============================================================================
# Pool step-0 background replicates -> CDR3-agnostic background bg(g)
#
# Reads the folded replicates from generate_step0_background.R and, per chain,
# POOLS counts across replicates:
#   bg(g) = sum_r n_pass(g, rep_r) / sum_r n_all(g, rep_r)
# (a count-weighted average of the per-replicate pass rates). Also reports the
# per-replicate spread so you can see how stable each estimate is.
#
# Writes step0_background/bg_<label>_<chain>.csv with columns:
#   level, n_all, n_pass, bg_rate, rate_min, rate_max, n_reps
# for level = each V gene, each J gene, and each CDR3 length ("L_<n>").
#
# If the cognate step-0 folds exist (TCR_motif_atlas/<epitope>/step0), prints a
# side-by-side of the pooled decoy bg vs each cognate pass rate for TRAV12-2 and
# the CDR3-length gradient — the artifact check (epitope-independent => artifact).
# =============================================================================

.sourced_for_benchmark <- TRUE
source("ESM_motif_builder.R")

OUT_DIR     <- "step0_background"
BG_PEPTIDES <- c("A0201_ALAAAAAAV")
COGNATE     <- epitopes                 # from config: GIL, LLW
THRESHOLD   <- ESM_THRESHOLDS[1]        # 0.5 = the step-1 selection threshold
FP_V        <- "TRAV12-2"               # false-positive V to spotlight (alpha)

# per-(level) n_all / n_pass for one folded step-0 dir + chain, keyed by V, J, L
level_counts <- function(step0_dir, chain_letter, threshold = THRESHOLD) {
  chain_name  <- if (chain_letter == "A") "alpha" else "beta"
  scores_file <- file.path(step0_dir, sprintf("output_%s.csv", chain_name))
  model_file  <- file.path(step0_dir, sprintf("model_%s.csv",  chain_name))
  if (!file.exists(scores_file) || !file.exists(model_file)) return(NULL)
  s0   <- load_esm_scores(scores_file, model_file)
  pass <- s0[[SCORE_COL]] >= threshold

  one <- function(key, type) {
    keep <- !is.na(key)
    a <- table(key[keep]); p <- table(key[keep & pass])
    lv <- names(a)
    data.frame(type = type,
               level = lv,
               n_all  = as.integer(a[lv]),
               n_pass = as.integer(ifelse(lv %in% names(p), p[lv], 0)),
               stringsAsFactors = FALSE)
  }
  rbind(
    one(s0[[paste0("TR", chain_letter, "V")]], "V"),
    one(s0[[paste0("TR", chain_letter, "J")]], "J"),
    one(paste0("L_", nchar(s0[[paste0("cdr3_TR", chain_letter)]])), "length")
  )
}

# pool a level-count table across replicate dirs for one label + chain
pool_background <- function(label, chain_letter) {
  reps <- list.dirs(OUT_DIR, recursive = FALSE)
  reps <- reps[grepl("rep[0-9]+$", basename(reps))]
  per_rep <- Filter(Negate(is.null),
                    lapply(reps, function(rd) level_counts(file.path(rd, label, "step0"), chain_letter)))
  if (length(per_rep) == 0) return(NULL)

  all <- do.call(rbind, Map(function(df, i) { df$rep <- i; df }, per_rep, seq_along(per_rep)))
  all$rate <- all$n_pass / all$n_all

  agg <- aggregate(cbind(n_all, n_pass) ~ type + level, all, sum)
  agg$bg_rate <- round(agg$n_pass / agg$n_all, 4)
  spread <- aggregate(rate ~ type + level, all, function(x) c(min = min(x), max = max(x), n = length(x)))
  spread <- data.frame(type = spread$type, level = spread$level,
                       rate_min = round(spread$rate[, "min"], 4),
                       rate_max = round(spread$rate[, "max"], 4),
                       n_reps   = as.integer(spread$rate[, "n"]),
                       stringsAsFactors = FALSE)
  merge(agg, spread, by = c("type", "level"))
}

for (label in BG_PEPTIDES) {
  for (ch in c("A", "B")) {
    bg <- pool_background(label, ch)
    chain_name <- if (ch == "A") "alpha" else "beta"
    if (is.null(bg)) {
      message(sprintf("[%s | %s] no folded replicates yet.", label, chain_name)); next
    }
    # separate CSVs per level type: V genes, J genes, CDR3 lengths
    for (ty in c("V", "J", "length")) {
      sub <- bg[bg$type == ty, setdiff(names(bg), "type"), drop = FALSE]
      if (nrow(sub) == 0) next
      if (ty == "length") {
        sub$len <- as.integer(sub("^L_", "", sub$level))
        sub <- sub[order(sub$len), ]
      } else {
        sub <- sub[order(-sub$bg_rate), ]     # genes: strongest artifact first
      }
      out <- file.path(OUT_DIR, sprintf("bg_%s_%s_%s.csv", label, chain_name, ty))
      write.csv(sub, out, row.names = FALSE)
    }
    message(sprintf("[%s | %s] pooled %d replicate(s) -> bg_%s_%s_{V,J,length}.csv",
                    label, chain_name, max(bg$n_reps), label, chain_name))
  }
}

# ---- artifact check: pooled decoy bg vs cognate pass rates (alpha) ----------
cognate_rate <- function(epitope, chain_letter, level_prefix_fun) {
  lc <- level_counts(file.path("TCR_motif_atlas", epitope, "step0"), chain_letter)
  if (is.null(lc)) return(NULL)
  lc$rate <- round(lc$n_pass / lc$n_all, 3)
  lc
}

bgA <- pool_background(BG_PEPTIDES[1], "A")
if (!is.null(bgA)) {
  cat(sprintf("\n===== Artifact check (alpha, threshold %.2f): decoy bg vs cognate =====\n", THRESHOLD))

  cat(sprintf("\n%s pass rate:\n", FP_V))
  row <- data.frame(source = paste0("decoy(", BG_PEPTIDES[1], ")"),
                    rate = bgA$bg_rate[bgA$level == FP_V])
  for (ep in COGNATE) {
    cr <- cognate_rate(ep, "A")
    if (!is.null(cr)) row <- rbind(row, data.frame(source = ep, rate = cr$rate[cr$level == FP_V]))
  }
  print(row, row.names = FALSE)

  cat("\nCDR3-alpha length gradient (pass rate):\n")
  lens <- paste0("L_", 10:16)
  tab  <- data.frame(len = 10:16, decoy = bgA$bg_rate[match(lens, bgA$level)])
  for (ep in COGNATE) {
    cr <- cognate_rate(ep, "A")
    if (!is.null(cr)) tab[[ep]] <- cr$rate[match(lens, cr$level)]
  }
  print(tab, row.names = FALSE)
  cat("\n(If decoy ~ cognate for a level, that level's pass rate is epitope-independent => fold artifact.)\n")
}
