# =============================================================================
# Does continuous evidence (SOFT_VJ_EVIDENCE) give a better V/J posterior?
#
# Ground truth needs no folding: each epitope has a validation.csv whose
# Label == 1 rows are real, experimentally determined binders. A good posterior
# should put mass on the V and J genes those binders actually use.
#
# Metric: posterior mass on the true genes, where "true" is weighted by how
# often the positives use each gene. This is the probability that a gene drawn
# from the posterior matches a gene drawn from the real binders, i.e.
#     overlap = sum_g  posterior(g) * freq_positives(g)
# Higher is better. Reported as a ratio to the repertoire prior, so 1.0 means
# "no better than drawing genes by baseline frequency".
#
# Also reports the rank of the positives' single most used gene.
#
# Usage:  Rscript test_soft_vj.R
# =============================================================================

.sourced_for_benchmark <- TRUE
source("ESM_motif_builder.R")

RUN      <- "TCR_motif_atlas_cdr3_templated_new_installation_reduced_cycles"
EPITOPES <- c("A0201_ELAGIGILTV", "A0201_GILGFVFTL", "A0201_LLWNGPMAV")
STEPS    <- c(0, 1)                       # posterior is built FROM these steps
THRESH   <- ESM_THRESHOLDS                # step s uses ESM_THRESHOLDS[s+1]
T_GRID   <- c(0, 0.25, 0.5, 1.0)          # 0 = hard threshold (original)

overlap <- function(post, truth_freq) sum(post * truth_freq[names(post)], na.rm = TRUE)

rows <- list()
for (ep in EPITOPES) {
  val <- read.csv(file.path(RUN, ep, "validation.csv"), stringsAsFactors = FALSE)
  pos <- val[val$Label == 1, ]

  for (st in STEPS) {
    d <- file.path(RUN, ep, sprintf("step%d", st))
    for (ch in c("A", "B")) {
      cn <- if (ch == "A") "alpha" else "beta"
      sf <- file.path(d, sprintf("output_%s.csv", cn))
      mf <- meta_csv(d, cn)
      if (!file.exists(sf) || !file.exists(mf)) next
      scored <- load_esm_scores(sf, mf)
      thr    <- THRESH[st + 1]
      top    <- scored[scored[[SCORE_COL]] >= thr, ]
      if (nrow(top) < 10) next

      for (gt in c("V", "J")) {
        col <- paste0("TR", ch, gt)
        if (!col %in% colnames(pos)) next
        tf <- table(pos[[col]][!is.na(pos[[col]])])
        tf <- tf / sum(tf)                                   # positives' gene usage
        prior <- vj_baseline_prior[[ch]][[gt]]
        base  <- overlap(prior / sum(prior), tf)             # repertoire-prior reference

        for (tsd in T_GRID) {
          SOFT_VJ_EVIDENCE <<- tsd > 0
          SOFT_VJ_T_SD     <<- tsd
          p <- suppressMessages(
                 extract_vj_marginal_posterior(top, scored, ch,
                                               vj_baseline_prior[[ch]],
                                               threshold = thr))
          if (is.null(p)) next
          v <- p[[gt]]
          top_true <- names(tf)[which.max(tf)]
          rows[[length(rows) + 1L]] <- data.frame(
            epitope = sub("^A0201_", "", ep), step = st, gene = col,
            T_sd    = tsd,
            lift    = overlap(v, tf) / base,                 # vs repertoire prior
            rank_top_true = if (top_true %in% names(v))
                              which(order(-v) == which(names(v) == top_true)) else NA,
            n_genes = length(v), stringsAsFactors = FALSE)
        }
      }
    }
  }
}

res <- do.call(rbind, rows)
write.csv(res, "test_soft_vj_results.csv", row.names = FALSE)

cat("\n=== posterior mass on the true genes, relative to the repertoire prior ===\n")
cat("    (lift > 1 = better than baseline gene frequencies; higher is better)\n\n")
w <- reshape(res[, c("epitope","step","gene","T_sd","lift")],
             idvar = c("epitope","step","gene"), timevar = "T_sd", direction = "wide")
colnames(w) <- sub("^lift\\.", "T=", colnames(w))
print(w, row.names = FALSE, digits = 3)

cat("\n=== mean lift by T (T=0 is the current hard threshold) ===\n")
m <- aggregate(lift ~ T_sd, res, mean)
m$vs_hard <- round(m$lift / m$lift[m$T_sd == 0], 3)
print(m, row.names = FALSE, digits = 4)

cat("\n=== rank of the positives' most used gene ===\n")
r <- aggregate(rank_top_true ~ T_sd, res, function(x) mean(x, na.rm = TRUE))
print(r, row.names = FALSE, digits = 3)

cat("\n=== win/loss vs hard threshold, per epitope x step x gene ===\n")
for (tsd in setdiff(T_GRID, 0)) {
  a <- res[res$T_sd == tsd, ]; b <- res[res$T_sd == 0, ]
  k <- merge(a, b, by = c("epitope","step","gene"), suffixes = c("", "_hard"))
  cat(sprintf("  T=%.2f : better in %d of %d comparisons (median change %+.1f%%)\n",
              tsd, sum(k$lift > k$lift_hard), nrow(k),
              100 * (median(k$lift / k$lift_hard) - 1)))
}
