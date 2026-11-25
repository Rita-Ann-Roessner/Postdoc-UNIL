# =============================================================================
# Does scoring V/J COMBINATIONS beat scoring genes?
#
# Same static test as test_soft_vj.R, but the unit is now the V/J combination.
# Ground truth: the V/J combinations the validated positives actually use.
# Metric: sum_c P_sample(c) * freq_positives(c), divided by the same quantity
# under the repertoire prior. Lift > 1 means better than baseline frequencies.
#
# Compared:
#   prior       P(c) ∝ P_base(c)                                   (reference)
#   gene        P(c) ∝ post_V(v) * post_J(j)        <- CURRENT pipeline
#   cell        P(c) ∝ P_base(c) * exp(beta * z(s_hat(c)))         <- proposed
#     s_hat(c) = (n_c*mean_score(c) + tau*s_additive(c)) / (n_c + tau)
#     s_additive = fitted values of lm(score ~ V + J), the factorized model
#                  written in score space
#     tau        = within-combination var / between-combination var
#
# beta is swept, NOT fitted. All candidate combinations are those with a CDR3
# baseline, exactly as sample_vj_pairs restricts them.
#
# Usage:  Rscript test_cell_vj.R
# =============================================================================

.sourced_for_benchmark <- TRUE
source("ESM_motif_builder.R")

RUN      <- "TCR_motif_atlas_cdr3_templated_new_installation_reduced_cycles"
EPITOPES <- c("A0201_ELAGIGILTV", "A0201_GILGFVFTL", "A0201_LLWNGPMAV")
STEPS    <- c(0, 1)
BETAS    <- c(0.5, 1, 2, 4)
TAU_FIX  <- NA        # NA = estimate per chain/step from variance components

rows <- list()
for (ep in EPITOPES) {
  val <- read.csv(file.path(RUN, ep, "validation.csv"), stringsAsFactors = FALSE)
  pos <- val[val$Label == 1, ]

  for (st in STEPS) {
    d   <- file.path(RUN, ep, sprintf("step%d", st))
    thr <- ESM_THRESHOLDS[st + 1]

    for (ch in c("A", "B")) {
      cn <- if (ch == "A") "alpha" else "beta"
      sf <- file.path(d, sprintf("output_%s.csv", cn)); mf <- meta_csv(d, cn)
      if (!file.exists(sf) || !file.exists(mf)) next
      sc <- load_esm_scores(sf, mf)
      if (sum(sc[[SCORE_COL]] >= thr, na.rm = TRUE) < 10) next

      v_col <- paste0("TR", ch, "V"); j_col <- paste0("TR", ch, "J")
      sc$cell <- paste(sc[[v_col]], sc[[j_col]], sep = "_")
      y <- sc[[SCORE_COL]]

      # candidate combinations: those with a CDR3 baseline AND functional genes
      keys <- names(cdr3_baseline[[paste0("TR", ch)]])
      spl  <- strsplit(keys, "_", fixed = TRUE)
      vg   <- vapply(spl, `[`, character(1), 1); jg <- vapply(spl, `[`, character(1), 2)
      prior <- vj_baseline_prior[[ch]]
      ok   <- vg %in% names(prior$V) & jg %in% names(prior$J)
      keys <- keys[ok]; vg <- vg[ok]; jg <- jg[ok]

      # truth: positives' combination usage, restricted to candidates
      pc <- paste(pos[[v_col]], pos[[j_col]], sep = "_")
      pc <- pc[pc %in% keys]
      if (length(pc) < 10) next
      tf <- table(pc) / length(pc)

      p_base <- as.numeric(prior$V[vg]) * as.numeric(prior$J[jg])
      p_base <- p_base / sum(p_base)
      truth  <- as.numeric(tf[keys]); truth[is.na(truth)] <- 0
      base_lift <- sum(p_base * truth)

      # --- CURRENT: gene-level posterior ------------------------------------
      top <- sc[sc[[SCORE_COL]] >= thr, ]
      g   <- extract_vj_marginal_posterior(top, sc, ch, prior, threshold = thr)
      p_gene <- as.numeric(g$V[vg]) * as.numeric(g$J[jg]); p_gene <- p_gene / sum(p_gene)

      # --- PROPOSED: per-combination shrunk score ---------------------------
      n_c <- table(sc$cell); m_c <- tapply(y, sc$cell, mean)
      fit <- lm(y ~ sc[[v_col]] + sc[[j_col]])
      s_add_row <- fitted(fit)
      s_add_c   <- tapply(s_add_row, sc$cell, mean)

      if (is.na(TAU_FIX)) {
        rep_cells <- names(n_c)[n_c >= 3]
        if (length(rep_cells) >= 5) {
          vw <- mean(tapply(y, sc$cell, var)[rep_cells], na.rm = TRUE)
          vb <- max(var(m_c[rep_cells], na.rm = TRUE) - vw / mean(n_c[rep_cells]), 1e-6)
          tau <- vw / vb
        } else tau <- 2.2          # no repeats at this step: NY-ESO-1 estimate
      } else tau <- TAU_FIX

      nn <- as.numeric(n_c[keys]); nn[is.na(nn)] <- 0
      mm <- as.numeric(m_c[keys])
      aa <- as.numeric(s_add_c[keys])
      # combinations never folded: fall back to the additive prediction for the
      # pair, built from the gene-level fit (always defined)
      cf <- coef(fit); ic <- cf[1]
      gv <- cf[match(paste0("sc[[v_col]]", vg), names(cf))]; gv[is.na(gv)] <- 0
      gj <- cf[match(paste0("sc[[j_col]]", jg), names(cf))]; gj[is.na(gj)] <- 0
      aa[is.na(aa)] <- (ic + gv + gj)[is.na(aa)]
      mm[is.na(mm)] <- aa[is.na(mm)]

      s_hat <- (nn * mm + tau * aa) / (nn + tau)
      z <- (s_hat - mean(s_hat)) / stats::sd(s_hat)

      for (b in BETAS) {
        p_cell <- p_base * exp(b * z); p_cell <- p_cell / sum(p_cell)
        rows[[length(rows) + 1L]] <- data.frame(
          epitope = sub("^A0201_", "", ep), step = st, chain = cn,
          n_pos_cells = length(unique(pc)), tau = round(tau, 2), beta = b,
          gene = sum(p_gene * truth) / base_lift,
          cell = sum(p_cell * truth) / base_lift,
          stringsAsFactors = FALSE)
      }
    }
  }
}

res <- do.call(rbind, rows)
write.csv(res, "test_cell_vj_results.csv", row.names = FALSE)

cat("\n=== lift over the repertoire prior, per V/J COMBINATION ===\n")
cat("    gene = current pipeline, cell = proposed. Higher is better.\n\n")
print(res[res$beta == 1, c("epitope","step","chain","n_pos_cells","tau","gene","cell")],
      row.names = FALSE, digits = 3)

cat("\n=== mean over all comparisons, by beta ===\n")
a <- aggregate(cbind(gene, cell) ~ beta, res, mean)
a$ratio <- round(a$cell / a$gene, 3)
a$wins  <- aggregate(I(cell > gene) ~ beta, res, sum)[, 2]
a$n     <- aggregate(I(cell > gene) ~ beta, res, length)[, 2]
print(a, row.names = FALSE, digits = 4)
