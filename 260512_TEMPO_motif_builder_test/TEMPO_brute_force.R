library(TEMPOtrain)
library(pROC)

# =============================================================================
# Brute-force benchmark
#
# 1. Generate N_RANDOM random TCRs from the MixTCRviz baseline (once, shared
#    across all epitopes).
# 2. For each epitope:
#    a. Score all N_RANDOM TCRs with the cached TEMPO predictor.
#    b. Keep top TOP_PERC fraction (default 0.1% = 1,000 out of 1,000,000).
#    c. Train a TEMPO motif on those top TCRs.
#    d. Score validation.csv → AUC / AUC0.1.
# 3. Save results to auc_summary_brute.csv.
# =============================================================================

.sourced_by_optimizer <- TRUE
source("TEMPO_motif_builder_test.R", local = FALSE)

BRUTE_OUTPUT_DIR <- "IMMREP23/brute_force_runs_small" # "brute_force_runs"
N_RANDOM         <- 8267 #1e6    # random TCRs to generate
TOP_PERC         <- 0.001  # top fraction to keep (0.001 = 0.1%)

dir.create(BRUTE_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# Bulk random TCR generator (efficient: batch-draw CDR3s per unique V/J key)
# =============================================================================

generate_random_tcrs_bulk <- function(chain, n_total, cdr3_baseline) {
  chain_letter     <- if (chain == "TRA") "A" else "B"
  v_col            <- paste0("TR", chain_letter, "V")
  j_col            <- paste0("TR", chain_letter, "J")
  cdr3_col         <- paste0("cdr3_TR", chain_letter)
  valid_keys       <- names(cdr3_baseline[[chain]])

  # Sample V/J keys uniformly with replacement
  sampled_keys     <- sample(valid_keys, size = n_total, replace = TRUE)
  key_counts       <- table(sampled_keys)

  rows <- lapply(names(key_counts), function(key) {
    n_draw      <- as.integer(key_counts[[key]])
    vj_counts   <- cdr3_baseline[[chain]][[key]]
    lens_avail  <- names(vj_counts)

    # Length probabilities
    lens_cnt <- sapply(lens_avail, function(l) sum(vj_counts[[l]][, 1], na.rm = TRUE))
    if (sum(lens_cnt) == 0) return(NULL)

    sel_lens <- sample(lens_avail, size = n_draw, replace = TRUE,
                       prob = lens_cnt / sum(lens_cnt))

    seqs <- vapply(sel_lens, function(len) {
      count_mat <- as.matrix(vj_counts[[len]])
      prob_mat  <- apply(count_mat, 2, function(col) {
        if (sum(col) == 0) rep(1 / nrow(count_mat), nrow(count_mat))
        else               col / sum(col)
      })
      aa      <- rownames(prob_mat)
      seq_vec <- vapply(seq_len(ncol(prob_mat)), function(pos) {
        p <- prob_mat[, pos]
        if (sum(p) == 0) return(NA_character_)
        sample(aa, 1, prob = p)
      }, character(1))
      if (any(is.na(seq_vec))) return(NA_character_)
      paste0(seq_vec, collapse = "")
    }, character(1))

    vj_split <- strsplit(key, "_")[[1]]
    data.frame(
      V    = vj_split[1],
      J    = vj_split[2],
      cdr3 = seqs,
      stringsAsFactors = FALSE
    )
  })

  df <- do.call(rbind, rows)
  df <- df[!is.na(df$cdr3), ]
  colnames(df) <- c(v_col, j_col, cdr3_col)
  df
}


# =============================================================================
# Step 1: Generate random TCRs (once, shared across all epitopes)
# =============================================================================

random_csv <- file.path(BRUTE_OUTPUT_DIR, "random.csv")

if (file.exists(random_csv)) {
  message(sprintf("Loading cached random TCR pool from: %s", random_csv))
  random_tcrs <- read.csv(random_csv)
} else {
  message(sprintf("Generating %s random TCRs from MixTCRviz baseline...", format(N_RANDOM, big.mark = ",")))

  df_A <- generate_random_tcrs_bulk("TRA", N_RANDOM, cdr3_baseline)
  df_B <- generate_random_tcrs_bulk("TRB", N_RANDOM, cdr3_baseline)

  message(sprintf("  Alpha chain: %s TCRs", format(nrow(df_A), big.mark = ",")))
  message(sprintf("  Beta chain:  %s TCRs", format(nrow(df_B), big.mark = ",")))

  # Random pairing — align to shorter chain
  n_pair <- min(nrow(df_A), nrow(df_B))
  df_A   <- df_A[sample(nrow(df_A), n_pair), ]
  df_B   <- df_B[sample(nrow(df_B), n_pair), ]
  rownames(df_A) <- rownames(df_B) <- NULL

  random_tcrs        <- cbind(df_A, df_B)
  random_tcrs$id     <- seq_len(nrow(random_tcrs))

  write.csv(random_tcrs, random_csv, row.names = FALSE)
  message(sprintf("  Saved %s paired TCRs to: %s",
                  format(nrow(random_tcrs), big.mark = ","), random_csv))
}


# =============================================================================
# Step 2: Per-epitope scoring + top selection + validation
# =============================================================================

auc_summary <- data.frame()

for (epitope in epitopes) {
  mhc     <- sub("_.*", "", epitope)
  peptide <- sub("^[^_]*_", "", epitope)
  message(sprintf("\n===== %s =====", epitope))

  epitope_dir   <- file.path(BRUTE_OUTPUT_DIR, epitope)
  dir.create(epitope_dir, showWarnings = FALSE, recursive = TRUE)

  predictor_rds   <- file.path(BASE_OUTPUT_DIR, epitope, paste0("predictor_", epitope, ".rds"))
  validation_file <- file.path(BASE_OUTPUT_DIR, epitope, "validation.csv")

  if (!file.exists(predictor_rds)) {
    message(sprintf("  Predictor not found — building from %s/%s.csv", epitope, epitope))
    predictor_rds <- build_tempo_predictor(
      file.path(BASE_OUTPUT_DIR, epitope, paste0(epitope, ".csv")), predictor_rds
    )
  }

  tryCatch({
    # (a) Add model/peptide/MHC columns required by TEMPOtrain
    model_csv <- file.path(epitope_dir, "random_model.csv")
    if (!file.exists(model_csv)) {
      df_model           <- random_tcrs
      df_model$model     <- paste0(mhc, "_", peptide)
      df_model$peptide   <- peptide
      df_model$MHC       <- mhc
      write.csv(df_model, model_csv, row.names = FALSE)
    }

    # (b) Score all random TCRs
    message("  Scoring random TCR pool...")
    score_dir  <- file.path(epitope_dir, "scoring")
    pred_file  <- run_tempo_scoring(predictor_rds, model_csv, score_dir)
    scored     <- read.csv(pred_file)

    # (c) Keep top TOP_PERC fraction by perc_rank (lower = better)
    n_top      <- max(1L, round(nrow(scored) * TOP_PERC))
    top_scored <- scored[order(scored[[SCORE_COL]]), ][seq_len(n_top), ]
    message(sprintf("  Keeping top %d / %s TCRs (%.1f%%)",
                    n_top, format(nrow(scored), big.mark = ","), TOP_PERC * 100))

    top_csv <- file.path(epitope_dir, "top_tcrs.csv")
    write.csv(top_scored, top_csv, row.names = FALSE)

    # (d) Train TEMPO motif on top TCRs, score validation set
    message("  Validating motif on validation set...")
    val_result <- validate_motif(
      motif_file      = top_csv,
      validation_file = validation_file,
      output_dir      = file.path(epitope_dir, "validation"),
      peptide         = peptide,
      mhc             = mhc,
      score_col       = "score",
      plot            = FALSE
    )

    auc_summary <- rbind(auc_summary, data.frame(
      epitope  = epitope,
      n_top    = n_top,
      auc      = val_result$auc,
      auc01    = val_result$auc01,
      stringsAsFactors = FALSE
    ))

    write.csv(auc_summary,
              file.path(BRUTE_OUTPUT_DIR, "auc_summary_brute.csv"),
              row.names = FALSE)

  }, error = function(e) {
    message(sprintf("  ERROR: %s", e$message))
  })
}


message("\n===== Brute-force AUC Summary =====")
for (i in seq_len(nrow(auc_summary))) {
  message(sprintf("  %s: AUC = %.4f  |  AUC0.1 = %.4f",
                  auc_summary$epitope[i], auc_summary$auc[i], auc_summary$auc01[i]))
}
message("  Saved to: auc_summary_brute.csv")
