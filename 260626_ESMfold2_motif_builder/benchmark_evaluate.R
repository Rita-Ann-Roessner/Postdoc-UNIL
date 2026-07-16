# =============================================================================
# Benchmark evaluator — walks the grid trees produced by benchmark_generate.R
# (and folded on the cluster) and, for each (grid point, step), computes:
#   - n_top per chain
#   - motif AUC / AUC0.1 (TEMPO-validated against validation.csv)
#   - false-positive gene share in beta (default TRBV25-1) and the dominant
#     beta V gene + its share  -> the trajectory that tells you where the prior
#     stops the false positive from growing across steps.
# Writes benchmark/summary.csv and prints it.
# =============================================================================

.sourced_for_benchmark <- TRUE
source("ESM_motif_builder.R")   # functions + config (ESM_THRESHOLDS, epitopes, SCORE_COL)

BENCH_DIR  <- "benchmark"
STEPS      <- c(1,2)           # steps to evaluate (those you have folded)
FP_GENE    <- "TRBV25-1"       # false-positive V gene to track (beta chain)
PLOT_MOTIF <- TRUE             # also write a MixTCRviz motif logo per grid point / step / chain

# Step k's top binders are defined at ESM_THRESHOLDS[k+1] (the threshold the
# pipeline applies to step k's scores when seeding step k+1), so evaluating up
# to step max(STEPS) needs index max(STEPS)+1 in the schedule.
if (max(STEPS) + 1 > length(ESM_THRESHOLDS))
  stop(sprintf("ESM_THRESHOLDS has %d entries, but evaluating up to step %d needs index %d ",
               length(ESM_THRESHOLDS), max(STEPS), max(STEPS) + 1),
       "(threshold for step k is ESM_THRESHOLDS[k+1]). Extend ESM_THRESHOLDS.")

top_at <- function(out_csv, model_csv, thr) {
  if (!file.exists(out_csv) || !file.exists(model_csv)) return(NULL)
  sc <- load_esm_scores(out_csv, model_csv)
  sc[sc[[SCORE_COL]] >= thr, ]
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
for (gp_dir in list.dirs(BENCH_DIR, recursive = FALSE)) {
  gp <- basename(gp_dir)
  if (!grepl("^vj[0-9]+_len[0-9]+$", gp)) next
  vj  <- as.integer(sub("^vj0*([0-9]+)_.*",  "\\1", gp))
  len <- as.integer(sub(".*_len0*([0-9]+)$", "\\1", gp))

  for (ep in epitopes) {
    mhc <- sub("_.*", "", ep); peptide <- sub("^[^_]*_", "", ep)
    val <- file.path(gp_dir, ep, "validation.csv")

    for (s in STEPS) {
      sd  <- file.path(gp_dir, ep, sprintf("step%d", s))
      thr <- ESM_THRESHOLDS[s + 1]   # step s's top binders = step-(s+1) threshold
      top_a <- top_at(file.path(sd, "output_alpha.csv"), file.path(sd, "model_alpha.csv"), thr)
      top_b <- top_at(file.path(sd, "output_beta.csv"),  file.path(sd, "model_beta.csv"),  thr)
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

      # motif logos (both chains): alpha batch -> real alpha (+ dummy beta),
      # beta batch -> real beta (+ dummy alpha). Title carries the VJ value;
      # PDFs land in <grid>/<epitope>/step<s>/motif_{alpha,beta}/.
      if (PLOT_MOTIF) {
        lab <- sprintf("%s VJ=%d LEN=%d", ep, vj, len)
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
        vj = vj, len = len, epitope = ep, step = s,
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
  message("No folded steps found under ", BENCH_DIR,
          " — run benchmark_generate.R, fold, gather, then re-run.")
} else {
  summary <- do.call(rbind, rows)
  summary <- summary[order(summary$epitope, summary$step, summary$vj), ]
  write.csv(summary, file.path(BENCH_DIR, "summary.csv"), row.names = FALSE)
  cat("\n===== Benchmark summary (VJ sweep) =====\n"); print(summary, row.names = FALSE)
  cat(sprintf("\nSaved to %s\n", file.path(BENCH_DIR, "summary.csv")))
}
