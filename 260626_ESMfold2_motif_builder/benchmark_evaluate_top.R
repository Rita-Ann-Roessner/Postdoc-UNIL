# =============================================================================
# Benchmark evaluator (top-N variant) — same as benchmark_evaluate.R, but defines
# top binders as the TOP_N highest-scoring TCRs per chain (fixed sample size)
# instead of an iptm_pair_mean threshold. Fixed n removes the sample-size
# confound between grid points (every grid point's motif is built from the same
# number of binders), so AUC / gene-shares are compared on equal footing.
#
# For each (grid point, step) it computes:
#   - n_top per chain (== TOP_N, or fewer if the step has < TOP_N folded TCRs)
#   - motif AUC / AUC0.1 (TEMPO-validated against validation.csv)
#   - false-positive gene share in beta (default TRBV25-1) + dominant beta V gene
#   - (optional) MixTCRviz motif logo per chain: <sd>/motif_{alpha,beta}.pdf
# Writes <EVAL_DIR>/summary_top<N>.csv and prints it.
# =============================================================================

.sourced_for_benchmark <- TRUE
source("ESM_motif_builder.R")   # functions + config (epitopes, SCORE_COL); main loop skipped

EVAL_DIR   <- "test"          # directory holding the grid trees (set to "benchmark" if that's where yours live)
STEPS      <- c(1, 2)         # steps to evaluate (those you have folded)
TOP_N      <- 50              # take the TOP_N highest-scoring TCRs per chain (threshold-free)
FP_GENE    <- "TRBV25-1"      # false-positive V gene to track (beta chain)
PLOT_MOTIF <- TRUE            # also write a MixTCRviz motif logo per grid point / step / chain

# Top-N highest scorers per chain (no threshold).
top_n <- function(out_csv, model_csv, n) {
  if (!file.exists(out_csv) || !file.exists(model_csv)) return(NULL)
  sc <- load_esm_scores(out_csv, model_csv)
  sc <- sc[order(sc[[SCORE_COL]], decreasing = TRUE), ]
  head(sc, n)
}

# plot_step_motif (via MixTCRviz) writes <sd>/motif_<chain>/<Model_*>.pdf; flatten
# that to a single <sd>/motif_<chain>.pdf and drop the subdir.
plot_motif_flat <- function(top, chain_letter, lab, s, sd) {
  chain_name <- if (chain_letter == "A") "alpha" else "beta"
  plot_step_motif(top, chain_letter, lab, s, sd)          # -> sd/motif_<chain>/*.pdf
  d    <- file.path(sd, sprintf("motif_%s", chain_name))
  pdfs <- list.files(d, pattern = "\\.pdf$", full.names = TRUE)
  if (length(pdfs) >= 1) {
    file.rename(pdfs[1], file.path(sd, sprintf("motif_%s.pdf", chain_name)))
    unlink(d, recursive = TRUE)
  }
}

rows <- list()
for (gp_dir in list.dirs(EVAL_DIR, recursive = FALSE)) {
  gp <- basename(gp_dir)
  if (!grepl("^vj[0-9]+_len[0-9]+$", gp)) next
  vj  <- as.integer(sub("^vj", "",     sub("_.*", "", gp)))   # "vj030_len020" -> 30
  len <- as.integer(sub(".*_len", "", gp))                    #                -> 20

  for (ep in epitopes) {
    mhc <- sub("_.*", "", ep); peptide <- sub("^[^_]*_", "", ep)
    val <- file.path(gp_dir, ep, "validation.csv")

    for (s in STEPS) {
      sd <- file.path(gp_dir, ep, sprintf("step%d", s))
      top_a <- top_n(file.path(sd, "output_alpha.csv"), file.path(sd, "model_alpha.csv"), TOP_N)
      top_b <- top_n(file.path(sd, "output_beta.csv"),  file.path(sd, "model_beta.csv"),  TOP_N)
      if (is.null(top_a) && is.null(top_b)) next   # step not folded yet

      # motif AUC (TEMPO-validated)
      auc <- NA_real_; auc01 <- NA_real_
      motif <- build_motif_df(top_a, top_b, paste0(mhc, "_", peptide))
      if (!is.null(motif) && file.exists(val)) {
        mfile <- file.path(sd, "motif_eval.csv"); write.csv(motif, mfile, row.names = FALSE)
        r <- tryCatch(
          validate_motif(mfile, val, file.path(sd, "eval"), peptide, mhc, plot = FALSE),
          error = function(e) { message("  [", gp, "|", ep, "|step", s, "] val err: ",
                                        conditionMessage(e)); NULL })
        if (!is.null(r)) { auc <- r$auc; auc01 <- r$auc01 }
      }

      # motif logos (both chains)
      if (PLOT_MOTIF) {
        lab <- sprintf("%s VJ=%d LEN=%d topN=%d", ep, vj, len, TOP_N)
        if (!is.null(top_a) && nrow(top_a) >= 10) plot_motif_flat(top_a, "A", lab, s, sd)
        if (!is.null(top_b) && nrow(top_b) >= 10) plot_motif_flat(top_b, "B", lab, s, sd)
      }

      # false-positive gene diagnostics (beta)
      fp_share <- NA_real_; topV <- NA_character_; topV_share <- NA_real_
      if (!is.null(top_b) && nrow(top_b) > 0) {
        fp_share   <- mean(top_b$TRBV == FP_GENE, na.rm = TRUE)
        tb         <- sort(table(top_b$TRBV), decreasing = TRUE)
        topV       <- names(tb)[1]
        topV_share <- as.numeric(tb[1]) / nrow(top_b)
      }

      rows[[length(rows) + 1]] <- data.frame(
        vj = vj, len = len, epitope = ep, step = s, top_n = TOP_N,
        n_top_a = if (is.null(top_a)) NA_integer_ else nrow(top_a),
        n_top_b = if (is.null(top_b)) NA_integer_ else nrow(top_b),
        auc = auc, auc01 = auc01,
        fp_gene = FP_GENE, fp_share_beta = round(fp_share, 3),
        top_beta_V = topV, top_beta_V_share = round(topV_share, 3),
        stringsAsFactors = FALSE
      )
    }
  }
}

if (length(rows) == 0) {
  message("No folded steps found under ", EVAL_DIR,
          " — check EVAL_DIR / that the grid trees contain output_<chain>.csv.")
} else {
  summary <- do.call(rbind, rows)
  summary <- summary[order(summary$epitope, summary$step, summary$vj), ]
  out <- file.path(EVAL_DIR, sprintf("summary_top%d.csv", TOP_N))
  write.csv(summary, out, row.names = FALSE)
  cat(sprintf("\n===== Benchmark summary (top-%d filtering) =====\n", TOP_N))
  print(summary, row.names = FALSE)
  cat(sprintf("\nSaved to %s\n", out))
}
