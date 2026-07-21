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

BENCH_DIR  <- "benchmark_MUT"
STEPS      <- c(1, 2)          # integer folded steps to evaluate
FILTER     <- "threshold"    # "threshold": every step uses ESM_THRESHOLDS[s+1] (as before).
                              # "percentile": the LAST step uses the Nth-percentile cutoff of its
                              #   pooled alpha+beta scores; earlier steps still use the threshold.
PERCENTILE <- 95              # percentile for the last step when FILTER == "percentile"
FP_GENE    <- "TRBV25-1"      # false-positive V gene to track (beta chain)
PLOT_MOTIF <- TRUE            # also write a MixTCRviz motif logo per grid point / step / chain

# Threshold-based steps need ESM_THRESHOLDS[s+1]. When FILTER == "percentile" the
# LAST step uses a percentile cutoff instead, so it needs no schedule entry.
thr_steps <- if (FILTER == "percentile") STEPS[STEPS != max(STEPS)] else STEPS
if (length(thr_steps) > 0 && max(thr_steps) + 1 > length(ESM_THRESHOLDS))
  stop(sprintf("ESM_THRESHOLDS has %d entries, but threshold steps up to %d need index %d.",
               length(ESM_THRESHOLDS), max(thr_steps), max(thr_steps) + 1))

# Full scored TCRs for a chain (NULL if not folded yet).
load_scored <- function(out_csv, model_csv) {
  if (!file.exists(out_csv) || !file.exists(model_csv)) return(NULL)
  load_esm_scores(out_csv, model_csv)
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
  if (!grepl("^vj[0-9]+_len[0-9]+_mut[0-9]+$", gp)) next
  vj  <- as.integer(sub("^vj0*([0-9]+)_.*",       "\\1", gp))
  len <- as.integer(sub(".*_len0*([0-9]+)_.*",    "\\1", gp))
  mut <- as.integer(sub(".*_mut0*([0-9]+)$",      "\\1", gp)) / 1000   # basis points -> fraction

  for (ep in epitopes) {
    mhc <- sub("_.*", "", ep); peptide <- sub("^[^_]*_", "", ep)
    val <- file.path(gp_dir, ep, "validation.csv")

    for (s in STEPS) {
      sd <- file.path(gp_dir, ep, sprintf("step%d", s))
      a  <- load_scored(file.path(sd, "output_alpha.csv"), file.path(sd, "model_alpha.csv"))
      b  <- load_scored(file.path(sd, "output_beta.csv"),  file.path(sd, "model_beta.csv"))
      if (is.null(a) && is.null(b)) next   # step not folded yet

      # Cutoff: percentile of this step's pooled alpha+beta scores (last step when
      # FILTER == "percentile"), otherwise the fixed step-(s+1) threshold.
      if (FILTER == "percentile" && s == max(STEPS)) {
        pooled <- c(if (!is.null(a)) a[[SCORE_COL]], if (!is.null(b)) b[[SCORE_COL]])
        thr    <- as.numeric(quantile(pooled, PERCENTILE / 100, na.rm = TRUE))
      } else {
        thr <- ESM_THRESHOLDS[s + 1]
      }
      top_a <- if (!is.null(a)) a[a[[SCORE_COL]] >= thr, ] else NULL
      top_b <- if (!is.null(b)) b[b[[SCORE_COL]] >= thr, ] else NULL

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
        lab <- sprintf("%s VJ=%d LEN=%d MUT=%g", ep, vj, len, mut)
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
        vj = vj, len = len, mut = mut, epitope = ep, step = s,
        cutoff = if (FILTER == "percentile" && s == max(STEPS)) sprintf("p%g", PERCENTILE) else "thr",
        thr = round(thr, 3),
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
  summary <- summary[order(summary$epitope, summary$step, summary$vj, summary$len, summary$mut), ]
  write.csv(summary, file.path(BENCH_DIR, "summary.csv"), row.names = FALSE)
  cat("\n===== Benchmark summary (VJ sweep) =====\n"); print(summary, row.names = FALSE)
  cat(sprintf("\nSaved to %s\n", file.path(BENCH_DIR, "summary.csv")))
}
