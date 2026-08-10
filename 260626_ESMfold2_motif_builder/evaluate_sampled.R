# =============================================================================
# Evaluate a MULTI-SAMPLE ESMFold re-fold of top binders.
#
# The <SAMPLE_SUB> dir (e.g. step2_sample) holds output_<chain>.csv with SEVERAL
# folds per TCR (multiple seeds / sample_idx) — so each ID appears many times.
# This script:
#   1. averages iptm_pair_mean per ID  (one robust score per TCR),
#   2. joins the per-ID score to the gene metadata (TRAV/TRAJ/cdr3...) from the
#      parent <META_SUB>/model_<chain>.csv (keyed on the same IDs),
#   3. applies an ESM threshold -> top binders, writes top_tcrs_<chain>.csv,
#   4. builds the combined motif (build_motif_df) and writes it,
#   5. validates the motif against validation.csv (TEMPO) -> AUC / AUC0.1.
#
# Reuses the pipeline functions; only the loading (averaging + metadata join)
# is new. Non-destructive: writes top_tcrs_<chain>.csv / motif_top_binders.csv /
# validation/ into the sample dir and does NOT touch the *_seqs.csv fold inputs.
# =============================================================================

.sourced_for_benchmark <- TRUE
source("ESM_motif_builder.R")   # build_motif_df, validate_motif, write_fold_input_csv, setup

# ---- configuration ----------------------------------------------------------
BASE_DIR   <- "TCR_motif_atlas_no_decoy_correction"
EPITOPES   <- c("A0201_ELAGIGILTV", "A0201_GILGFVFTL", "A0201_LLWNGPMAV")     # add more labels as needed
SAMPLE_SUB <- "step2_sample"           # dir with the multi-fold output_<chain>.csv
META_SUB   <- "step2"                  # dir with model_<chain>.csv (gene metadata)
THRESHOLD  <- 0.65                      # ESM cutoff on the PER-ID AVERAGED iptm_pair_mean
WRITE_SEQS <- FALSE                    # also write top_tcrs_<chain>_seqs.csv (fold.py format)?
PLOT       <- TRUE

# ---- averaged-score loader --------------------------------------------------
# Average score_col over the repeated folds per ID, then join to the gene
# metadata from the parent step's model table (inner join on id).
load_sampled_scores <- function(sample_out, meta_model, score_col = SCORE_COL) {
  out  <- read.csv(sample_out, check.names = FALSE, stringsAsFactors = FALSE)
  meta <- read.csv(meta_model, check.names = FALSE, stringsAsFactors = FALSE)
  fix_id <- function(df) { h <- which(tolower(names(df)) == "id"); if (length(h) >= 1) names(df)[h[1]] <- "id"; df }
  out <- fix_id(out); meta <- fix_id(meta)
  if (!score_col %in% names(out)) stop("Column '", score_col, "' not in ", sample_out)

  avg <- aggregate(out[[score_col]], by = list(id = out$id), FUN = mean)
  names(avg)[2] <- score_col
  nf  <- aggregate(out[[score_col]], by = list(id = out$id), FUN = length)
  names(nf)[2] <- "n_folds"
  merge(meta, merge(avg, nf, by = "id"), by = "id")   # keep sampled IDs that have metadata
}

# ---- per epitope ------------------------------------------------------------
auc_rows <- list()
for (ep in EPITOPES) {
  mhc     <- sub("_.*", "", ep)
  peptide <- sub("^[^_]*_", "", ep)
  s_dir   <- file.path(BASE_DIR, ep, SAMPLE_SUB)
  m_dir   <- file.path(BASE_DIR, ep, META_SUB)
  val     <- file.path(BASE_DIR, ep, "validation.csv")

  message(sprintf("\n========== %s | %s (avg over folds, threshold %.2f) ==========", ep, SAMPLE_SUB, THRESHOLD))

  collect_top <- function(chain_letter) {
    chain_name  <- if (chain_letter == "A") "alpha" else "beta"
    sample_out  <- file.path(s_dir, sprintf("output_%s.csv", chain_name))
    meta_model  <- file.path(m_dir, sprintf("model_%s.csv",  chain_name))
    if (!file.exists(sample_out)) { warning("missing ", sample_out); return(NULL) }
    if (!file.exists(meta_model)) { warning("missing ", meta_model); return(NULL) }

    scored <- load_sampled_scores(sample_out, meta_model)
    top    <- scored[scored[[SCORE_COL]] >= THRESHOLD, ]
    message(sprintf("  [chain %s] %d / %d IDs pass avg %s >= %.2f  (avg over ~%.1f folds/ID)",
                    chain_letter, nrow(top), nrow(scored), SCORE_COL, THRESHOLD, mean(scored$n_folds)))

    write.csv(top, file.path(s_dir, sprintf("top_tcrs_%s.csv", chain_name)), row.names = FALSE)
    if (WRITE_SEQS)
      write_fold_input_csv(top, file.path(s_dir, sprintf("top_tcrs_%s_avg_seqs.csv", chain_name)))
    top
  }

  top_alpha <- collect_top("A")
  top_beta  <- collect_top("B")

  motif_df <- build_motif_df(top_alpha, top_beta, paste0(mhc, "_", peptide))
  if (is.null(motif_df)) { warning(sprintf("[%s] no top binders — skipping motif/AUC.", ep)); next }
  motif_file <- file.path(s_dir, "motif_top_binders.csv")
  write.csv(motif_df, motif_file, row.names = FALSE)
  message(sprintf("  motif: %d TCRs -> %s", nrow(motif_df), motif_file))

  if (file.exists(val)) {
    res <- validate_motif(motif_file, val, file.path(s_dir, "validation"), peptide, mhc, plot = PLOT)
    message(sprintf("[%s] AUC = %.4f  |  AUC0.1 = %.4f", ep, res$auc, res$auc01))
    auc_rows[[ep]] <- data.frame(epitope = ep, threshold = THRESHOLD,
                                 n_alpha = if (is.null(top_alpha)) 0 else nrow(top_alpha),
                                 n_beta  = if (is.null(top_beta))  0 else nrow(top_beta),
                                 auc = res$auc, auc01 = res$auc01)
  } else {
    warning(sprintf("[%s] no validation.csv at %s — skipping AUC.", ep, val))
  }
}

if (length(auc_rows) > 0) {
  summary <- do.call(rbind, auc_rows)
  out_csv <- file.path(BASE_DIR, sprintf("sampled_auc_summary_%s.csv", SAMPLE_SUB))
  write.csv(summary, out_csv, row.names = FALSE)
  cat("\n===== Sampled AUC summary =====\n"); print(summary, row.names = FALSE)
  cat(sprintf("Saved to %s\n", out_csv))
}
