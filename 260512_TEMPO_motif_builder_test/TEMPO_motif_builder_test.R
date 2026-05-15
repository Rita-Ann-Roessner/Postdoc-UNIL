library(TEMPOtrain)
library(MixTCRviz)
library(dplyr)
library(tidyr)
library(pROC)

# =============================================================================
# TEMPO Motif Builder
#
# Three-stage refinement applied at each iteration:
#   (1) V/J enrichment     — correct the gene frame
#   (2) CDR3 length        — re-weight lengths from high-scorers
#   (3) PSSM blending      — bias CDR3 sequence content from high-scorers
#
# Iteration 0: random TCRs from the full V/J × CDR3 baseline (no priors)
# Iteration k: enriched V/J pairs, enriched CDR3 lengths, PSSM-blended CDR3s
# Final:       validate motif on labelled validation.csv → AUC
#
# NOTE: set SCORE_COL to the column name TEMPOtrain writes in its prediction
#       output (check head() of the pred file after the first run).
# =============================================================================

### ---- Configuration --------------------------------------------------------

INPUT_DIR       <- "/Users/roessner/Documents/PostDoc/Data/MixTCRviz/data_raw/HomoSapiens"
BASE_OUTPUT_DIR <- "."
SCORE_COL       <- "perc_rank" # lower = better binder (HLA convention)
N_ITERATIONS    <- 3           # scoring+enrichment steps (steps 1–N); step 0 is always flat random
N_ENRICHMENT    <- 19          # top-N V/J genes kept per column after enrichment
N_CDR3_MULTI    <- 5           # CDR3 sequences sampled per V/J pair (iterations > 0)
INIT_PERC_RANK  <- 10          # fixed threshold for scoring initial random pairs (pre-loop)
PERC_RANK_MID <- 9.8           # start of geometric decay schedule (optimizable; range 1–10)
                               # full schedule: [INIT_PERC_RANK, start, sqrt(start*0.05), 0.05]
PSSM_WEIGHT     <- 0.66        # blend weight: 0 = pure baseline, 1 = pure PSSM
MIN_TCRS_PSSM   <- 30         # minimum high-scorers required before using a PSSM

dico <- list(
  "LLWNGPMAV"  = "A0201",
  "ELAGIGILTV" = "A0201",
  "GILGFVFTL"  = "A0201"
)


# =============================================================================
# Module 1 — Iteration 0: random TCR generation from full V/J baseline
# =============================================================================

draw_random_cdr3 <- function(chain, v_seg, j_seg, cdr3_baseline) {
  key <- paste0(v_seg, "_", j_seg)
  if (!key %in% names(cdr3_baseline[[chain]])) return(NA)

  vj_counts         <- cdr3_baseline[[chain]][[key]]
  lengths_available <- names(vj_counts)
  lengths_cnt <- sapply(lengths_available, function(l) {
    sum(vj_counts[[l]][, 1], na.rm = TRUE)
  })
  if (sum(lengths_cnt) == 0) return(NA)

  len      <- sample(lengths_available, 1, prob = lengths_cnt / sum(lengths_cnt))
  prob_mat <- apply(as.matrix(vj_counts[[len]]), 2, function(col) {
    if (sum(col) == 0) rep(0, length(col)) else col / sum(col)
  })
  seq_vec <- sapply(seq_len(ncol(prob_mat)), function(pos) {
    p <- prob_mat[, pos]
    if (sum(p) == 0) return(NA)
    sample(rownames(prob_mat), 1, prob = p)
  })
  if (any(is.na(seq_vec))) return(NA)
  paste0(seq_vec, collapse = "")
}


generate_vj_pairs <- function(chain, input_dir, output_dir) {
  genes <- c(paste0("TR", chain, "V"), paste0("TR", chain, "J"))
  lst <- lapply(genes, function(gene) {
    df <- read.csv(file.path(input_dir, paste0(gene, ".csv")))
    df[df$gene_type == "F", 1]
  })
  pairs          <- expand.grid(lst[[1]], lst[[2]])
  colnames(pairs) <- genes
  write.csv(pairs, file.path(output_dir, paste0(genes[1], "_", genes[2], ".csv")),
            row.names = FALSE)
  pairs
}


sample_chain_cdr3 <- function(chain, pair_file, cdr3_baseline, output_file) {
  df           <- read.csv(pair_file)
  chain_letter <- sub("^TR", "", chain)   # "TRA" -> "A", "TRB" -> "B"
  v_col        <- paste0("TR", chain_letter, "V")
  j_col        <- paste0("TR", chain_letter, "J")
  df$random_CDR3 <- mapply(draw_random_cdr3,
                            chain = chain, v_seg = df[[v_col]], j_seg = df[[j_col]],
                            MoreArgs = list(cdr3_baseline = cdr3_baseline))
  df <- na.omit(df)
  write.csv(df, output_file, row.names = FALSE)
  df
}


pair_alpha_beta <- function(file_a, file_b, output_file, seed = 42) {
  df_A <- read.csv(file_a)
  df_B <- read.csv(file_b)
  colnames(df_A)[colnames(df_A) == "random_CDR3"] <- "cdr3_TRA"
  colnames(df_B)[colnames(df_B) == "random_CDR3"] <- "cdr3_TRB"

  n_A <- nrow(df_A); n_B <- nrow(df_B)
  if (n_A == 0 || n_B == 0) stop("Empty chain dataframe after CDR3 sampling.")
  n_max <- max(n_A, n_B)

  set.seed(seed)
  df_A <- df_A[sample(n_A, n_max, replace = n_A < n_max), , drop = FALSE]
  df_B <- df_B[sample(n_B, n_max, replace = n_B < n_max), , drop = FALSE]
  df_A <- df_A[sample(n_max), , drop = FALSE]
  df_B <- df_B[sample(n_max), , drop = FALSE]
  rownames(df_A) <- rownames(df_B) <- NULL

  df_paired <- cbind(df_A, df_B)
  write.csv(df_paired, output_file, row.names = FALSE)
  df_paired
}


prepare_model_csv <- function(df_paired, peptide, mhc, output_file) {
  df         <- df_paired
  df$id      <- seq_len(nrow(df))
  df$model   <- paste0(mhc, "_", peptide)
  df$peptide <- peptide
  df$MHC     <- mhc
  df$species <- "HomoSapiens"
  write.csv(df, output_file, row.names = FALSE)
  df
}


# =============================================================================
# Module 2 — PSSM building from high-scoring TCRs
#
# Pooled across all V/J pairs for a given chain (more data → more stable).
# The V/J-specific baseline already encodes V/J anchor patterns; the PSSM
# captures the variable positions that distinguish good binders.
# =============================================================================

build_cdr3_pssm <- function(cdr3_seqs, pseudocount = 0.1) {
  cdr3_seqs <- cdr3_seqs[!is.na(cdr3_seqs) & nchar(cdr3_seqs) > 0]
  if (length(cdr3_seqs) == 0) return(NULL)

  amino_acids <- c("A","C","D","E","F","G","H","I","K","L",
                   "M","N","P","Q","R","S","T","V","W","Y")

  by_length <- split(cdr3_seqs, nchar(cdr3_seqs))

  pssm_list <- lapply(names(by_length), function(len_str) {
    seqs <- by_length[[len_str]]
    L    <- as.integer(len_str)
    mat  <- matrix(pseudocount, nrow = length(amino_acids), ncol = L,
                   dimnames = list(amino_acids, NULL))
    for (seq in seqs) {
      chars <- strsplit(seq, "")[[1]]
      for (pos in seq_along(chars)) {
        aa <- chars[pos]
        if (aa %in% amino_acids) mat[aa, pos] <- mat[aa, pos] + 1
      }
    }
    apply(mat, 2, function(col) col / sum(col))
  })

  names(pssm_list) <- paste0("L_", names(by_length))
  pssm_list
}


# =============================================================================
# Module 3 — Multi-CDR3 sampling with all three enrichments
#
# draw_random_cdr3_multi handles all three refinements in one place:
#   - len_dist=NULL, cdr3_pssm=NULL → pure baseline  (iteration 0)
#   - len_dist given,  cdr3_pssm=NULL → (1)+(2)      (early iterations)
#   - len_dist given,  cdr3_pssm given → (1)+(2)+(3)  (later iterations)
# =============================================================================

draw_random_cdr3_multi <- function(chain, v_seg, j_seg, cdr3_baseline,
                                    n          = 5,
                                    len_dist   = NULL,
                                    cdr3_pssm  = NULL,
                                    pssm_weight = 0.5) {
  key <- paste0(v_seg, "_", j_seg)
  if (!key %in% names(cdr3_baseline[[chain]])) return(character(0))

  vj_counts         <- cdr3_baseline[[chain]][[key]]
  lengths_available <- names(vj_counts)

  # (2) Sample n lengths with replacement from enriched distribution
  if (!is.null(len_dist)) {
    valid_lens <- intersect(names(len_dist), lengths_available)
    if (length(valid_lens) > 0) {
      probs         <- unlist(len_dist[valid_lens])
      probs         <- probs / sum(probs)
      selected_lens <- sample(valid_lens, size = n, replace = TRUE, prob = probs)
    } else {
      # no overlap between enriched lengths and this V/J pair — sample uniformly
      selected_lens <- sample(lengths_available, size = min(n, length(lengths_available)))
    }
  } else {
    selected_lens <- sample(lengths_available, size = min(n, length(lengths_available)))
  }

  seqs <- vapply(selected_lens, function(len_name) {
    count_mat <- as.matrix(vj_counts[[len_name]])

    # V/J-baseline probability matrix
    prob_mat <- apply(count_mat, 2, function(col) {
      if (sum(col) == 0) rep(1 / nrow(count_mat), nrow(count_mat))
      else               col / sum(col)
    })

    # (3) Blend with PSSM when available and dimensions match
    if (!is.null(cdr3_pssm) && !is.null(cdr3_pssm[[len_name]])) {
      pssm      <- cdr3_pssm[[len_name]]
      aa_common <- intersect(rownames(prob_mat), rownames(pssm))
      if (length(aa_common) >= 10 && ncol(pssm) == ncol(prob_mat)) {
        prob_mat[aa_common, ] <- (1 - pssm_weight) * prob_mat[aa_common, ] +
                                       pssm_weight  * pssm[aa_common, ]
        prob_mat <- apply(prob_mat, 2, function(col) col / sum(col))
      }
    }

    amino_acids <- rownames(prob_mat)
    seq_vec <- sapply(seq_len(ncol(prob_mat)), function(pos) {
      p <- prob_mat[, pos]
      if (sum(p) == 0) return(NA_character_)
      sample(amino_acids, 1, prob = p)
    })
    if (any(is.na(seq_vec))) return(NA_character_)
    paste0(seq_vec, collapse = "")
  }, character(1))

  seqs[!is.na(seqs)]
}


sample_chain_cdr3_multi <- function(chain, pair_file, cdr3_baseline, output_file,
                                     n           = 5,
                                     len_dist    = NULL,
                                     cdr3_pssm   = NULL,
                                     pssm_weight = 0.5) {
  df           <- read.csv(pair_file)
  chain_letter <- sub("^TR", "", chain)   # "TRA" -> "A", "TRB" -> "B"
  v_col        <- paste0("TR", chain_letter, "V")
  j_col        <- paste0("TR", chain_letter, "J")

  results <- mapply(draw_random_cdr3_multi,
                    chain = chain, v_seg = df[[v_col]], j_seg = df[[j_col]],
                    MoreArgs = list(cdr3_baseline = cdr3_baseline,
                                   n           = n,
                                   len_dist    = len_dist,
                                   cdr3_pssm   = cdr3_pssm,
                                   pssm_weight = pssm_weight),
                    SIMPLIFY = FALSE)

  df_exp <- data.frame(
    V           = rep(df[[v_col]], lengths(results)),
    J           = rep(df[[j_col]], lengths(results)),
    random_CDR3 = unlist(results)
  )
  colnames(df_exp)[1:2] <- c(v_col, j_col)
  df_exp <- na.omit(df_exp)
  write.csv(df_exp, output_file, row.names = FALSE)
  df_exp
}


pair_alpha_beta_multi <- function(file_a, file_b, output_file, n = 5) {
  df_A <- read.csv(file_a)
  df_B <- read.csv(file_b)
  colnames(df_A)[colnames(df_A) == "random_CDR3"] <- "cdr3_TRA"
  colnames(df_B)[colnames(df_B) == "random_CDR3"] <- "cdr3_TRB"

  df_A <- df_A[rep(seq_len(nrow(df_A)), each = n), ]
  df_B <- df_B[rep(seq_len(nrow(df_B)), each = n), ]
  n_A  <- nrow(df_A); n_B <- nrow(df_B); n_max <- max(n_A, n_B)

  df_A <- df_A[sample(n_A, n_max, replace = n_A < n_max), ]
  df_B <- df_B[sample(n_B, n_max, replace = n_B < n_max), ]
  df_A <- df_A[sample(n_max), ]; df_B <- df_B[sample(n_max), ]
  rownames(df_A) <- rownames(df_B) <- NULL

  df_paired <- cbind(df_A, df_B)
  write.csv(df_paired, output_file, row.names = FALSE)
  df_paired
}


# =============================================================================
# Module 4a — MixTCRviz motif plot per iteration
# =============================================================================

plot_iter_motif <- function(top_tcrs, all_tcrs, peptide, iter, output_dir) {
  message(sprintf("  Plotting MixTCRviz motif (iter %d, n_top=%d)", iter, nrow(top_tcrs)))
  MixTCRviz(
    input1      = top_tcrs,
    input2      = all_tcrs,
    output.path = output_dir,
    plot        = TRUE,
    renormVJ    = TRUE,
    set.title   = sprintf("%s — iter %d top binders (n=%d)", peptide, iter, nrow(top_tcrs)),
    verbose     = 0
  )
}


# =============================================================================
# Module 4 — TEMPO predictor caching + scoring
#
# build_tempo_predictor():
#   Builds a TEMPO predictor from known binders once and saves result$predictor
#   as an RDS.  On subsequent calls the RDS is reused — skipping the slow
#   build.prank step.  The saved object is the raw predictor list (with $V, $J,
#   $L, $L_VJ, $CDR3, $CDR3_VJ) — NOT the full TEMPOtrain() return value —
#   because TEMPOtrain(input.train=<rds>) checks those fields directly.
#
# run_tempo_scoring():
#   Scores random TCRs using a cached predictor RDS (compute.prank=TRUE so the
#   output contains the "perc_rank" column used for iterative filtering).
# =============================================================================

build_tempo_predictor <- function(known_binders_file, predictor_rds) {
  if (file.exists(predictor_rds)) {
    message("  Predictor cache found — skipping rebuild: ", predictor_rds)
    return(predictor_rds)
  }
  message("  Building TEMPO predictor (build.prank=TRUE — this may take a few minutes)...")
  dir.create(dirname(predictor_rds), showWarnings = FALSE, recursive = TRUE)
  result <- TEMPOtrain(
    input.train      = known_binders_file,
    output.path      = dirname(predictor_rds),
    filename.train   = "predictor_train",
    build.prank      = TRUE,
    compute.prank    = FALSE,
    write.data.train = FALSE,
    write.predictor  = FALSE
  )
  saveRDS(result$predictor, predictor_rds)
  message("  Predictor saved to: ", predictor_rds)
  predictor_rds
}


run_tempo_scoring <- function(predictor_rds, pred_file, output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  TEMPOtrain(
    input.train      = predictor_rds,
    input.pred       = pred_file,
    output.path      = output_dir,
    filename.pred    = "pred",
    build.prank      = FALSE,
    compute.prank    = TRUE,
    write.data.pred  = TRUE,
    write.data.train = FALSE,
    write.predictor  = FALSE
  )
  pred_files <- list.files(output_dir, pattern = "^pred", full.names = TRUE)
  if (length(pred_files) == 0) stop("No prediction output found in ", output_dir)
  pred_files[1]
}


# =============================================================================
# Module 5 — Enrichment helpers
# =============================================================================

compute_vj_enrichment <- function(all_tcrs, top_tcrs,
                                   n    = NULL,
                                   cols = c("TRAV", "TRAJ", "TRBV", "TRBJ")) {
  count_df <- function(df, label) {
    df %>%
      select(all_of(cols)) %>%
      pivot_longer(everything(), names_to = "column", values_to = "value") %>%
      group_by(column, value) %>%
      summarise(count = n(), .groups = "drop") %>%
      rename(!!paste0(label, "_count") := count)
  }
  enrich <- count_df(all_tcrs, "baseline") %>%
    left_join(count_df(top_tcrs, "models"), by = c("column", "value")) %>%
    mutate(models_count = replace_na(models_count, 0),
           ratio        = models_count / baseline_count)

  if (!is.null(n)) {
    enrich <- enrich %>%
      group_by(column) %>%
      slice_max(order_by = ratio, n = n, with_ties = FALSE) %>%
      ungroup()
  }
  enrich
}


extract_cdr3_len_dist <- function(top_tcrs, chain) {
  cdr3_col <- paste0("cdr3_TR", chain)
  if (!cdr3_col %in% colnames(top_tcrs)) return(NULL)
  lens <- nchar(top_tcrs[[cdr3_col]])
  lens <- lens[!is.na(lens) & lens > 0]
  if (length(lens) == 0) return(NULL)
  tbl  <- table(paste0("L_", lens))
  as.list(tbl / sum(tbl))
}


load_tempo_scores <- function(pred_output_file, model_file, score_col) {
  scores <- read.csv(pred_output_file)
  model  <- read.csv(model_file)

  if (!score_col %in% colnames(scores))
    stop("Column '", score_col, "' not found in TEMPO output. ",
         "Available: ", paste(colnames(scores), collapse = ", "))

  if ("id" %in% colnames(scores) && "id" %in% colnames(model)) {
    merged <- merge(model, scores[, c("id", score_col)], by = "id")
  } else {
    if (nrow(scores) != nrow(model))
      stop("Row count mismatch between TEMPO output and model.csv; ",
           "ensure model.csv has a unique 'id' column.")
    merged        <- cbind(model, scores[[score_col]])
    colnames(merged)[ncol(merged)] <- score_col
  }
  merged
}


generate_enriched_vj_pairs <- function(enriched_df, output_dir) {
  for (chain in c("A", "B")) {
    genes <- c(paste0("TR", chain, "V"), paste0("TR", chain, "J"))
    lst   <- list(enriched_df$value[enriched_df$column == genes[1]],
                  enriched_df$value[enriched_df$column == genes[2]])
    pairs          <- expand.grid(lst[[1]], lst[[2]])
    colnames(pairs) <- genes
    write.csv(pairs, file.path(output_dir, paste0(genes[1], "_", genes[2], ".csv")),
              row.names = FALSE)
  }
}


# =============================================================================
# Module 6 — AUC validation
# Train TEMPO on the final motif; score validation.csv; compute AUC.
# =============================================================================

validate_motif <- function(motif_file, validation_file, output_dir,
                            peptide, mhc, score_col = "score") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  TEMPOtrain(
    input.train      = motif_file,
    input.pred       = validation_file,
    output.path      = output_dir,
    filename.train   = "motif_train",
    filename.pred    = "validation_pred",
    build.prank      = FALSE,
    compute.prank    = FALSE,
    write.data.pred  = TRUE,
    write.data.train = FALSE,
    write.predictor  = FALSE
  )

  pred_files <- list.files(output_dir, pattern = "^validation_pred", full.names = TRUE)
  if (length(pred_files) == 0) stop("No validation prediction output in ", output_dir)
  pred_df <- read.csv(pred_files[1])

  if (!"Label" %in% colnames(pred_df))
    stop("'Label' column not found in validation pred output; check that validation.csv contains a 'Label' column.")
  if (!score_col %in% colnames(pred_df))
    stop("Column '", score_col, "' not found in validation pred output.")

  pred_df <- pred_df[!is.na(pred_df[[score_col]]), ]

  roc_obj  <- roc(pred_df$Label, pred_df[[score_col]], quiet = TRUE)
  auc_val  <- as.numeric(auc(roc_obj))
  auc01_val <- as.numeric(auc(roc_obj,
                              partial.auc         = c(1, 0.9),
                              partial.auc.correct = TRUE,
                              partial.auc.focus   = "specificity"))

  plot(roc_obj,
       main = sprintf("ROC — %s_%s", mhc, peptide),
       col  = "#2166ac", lwd = 2,
       xaxt = "n",
       xlab = "FPR (1 - Specificity)",
       ylab = "TPR (Sensitivity)")
  axis(1, at = seq(1, 0, -0.2), labels = seq(0, 1, 0.2))
  legend("bottomright", bty = "o", bg = "white",
         box.col = "grey80",
         legend  = c(sprintf("AUC    = %.4f", auc_val),
                     sprintf("AUC0.1 = %.4f", auc01_val)))

  message(sprintf("  AUC = %.4f  |  AUC0.1 = %.4f", auc_val, auc01_val))
  list(auc = auc_val, auc01 = auc01_val, roc = roc_obj, pred = pred_df)
}


# =============================================================================
# Module 7 — Single enrichment step (score → filter → enrich → generate)
#
# Shared by the pre-loop step (threshold = INIT_PERC_RANK) and every loop
# iteration (threshold from geometric schedule).  Returns a list with the
# path to the newly generated model.csv and summary stats, or NULL if too
# few top binders were found.
# =============================================================================

run_enrichment_step <- function(predictor_rds, current_model_csv, threshold,
                                 step, peptide, mhc,
                                 base_output_dir, cdr3_baseline,
                                 n_enrichment, n_cdr3_multi,
                                 pssm_weight, min_tcrs_pssm,
                                 score_col, plot_motif = TRUE) {

  step_dir <- file.path(base_output_dir, sprintf("TEMPO_step%d_%s", step, peptide))
  dir.create(step_dir, showWarnings = FALSE, recursive = TRUE)

  # Score
  pred_output <- run_tempo_scoring(
    predictor_rds = predictor_rds,
    pred_file     = current_model_csv,
    output_dir    = file.path(step_dir, "TEMPO_scoring")
  )
  scored_tcrs <- load_tempo_scores(pred_output, current_model_csv, score_col)
  top_tcrs    <- scored_tcrs[scored_tcrs[[score_col]] < threshold, ]

  message(sprintf("  [Step %d] %d / %d TCRs with %s < %.2f",
                  step, nrow(top_tcrs), nrow(scored_tcrs), score_col, threshold))

  if (nrow(top_tcrs) < 10) {
    warning(sprintf("Too few high-scoring TCRs at step %d — stopping.", step))
    return(NULL)
  }

  if (plot_motif) {
    plot_iter_motif(top_tcrs, scored_tcrs, peptide, step,
                    output_dir = file.path(step_dir, "motif"))
  }

  # (1) V/J enrichment
  enriched   <- compute_vj_enrichment(scored_tcrs, top_tcrs, n = n_enrichment)

  # (2) CDR3 length distribution
  len_dist_A <- extract_cdr3_len_dist(top_tcrs, "A")
  len_dist_B <- extract_cdr3_len_dist(top_tcrs, "B")

  # (3) PSSM
  if (nrow(top_tcrs) >= min_tcrs_pssm) {
    pssm_A <- build_cdr3_pssm(top_tcrs$cdr3_TRA)
    pssm_B <- build_cdr3_pssm(top_tcrs$cdr3_TRB)
    message(sprintf("  PSSM built from %d high-scoring TCRs", nrow(top_tcrs)))
  } else {
    pssm_A <- pssm_B <- NULL
    message(sprintf("  Too few high-scorers for PSSM (%d < %d)", nrow(top_tcrs), min_tcrs_pssm))
  }

  # Generate enriched batch — saved in the same step_dir
  generate_enriched_vj_pairs(enriched, step_dir)

  sample_chain_cdr3_multi("TRA",
                          file.path(step_dir, "TRAV_TRAJ.csv"), cdr3_baseline,
                          file.path(step_dir, "TRAV_TRAJ_cdr3.csv"),
                          n = n_cdr3_multi, len_dist = len_dist_A,
                          cdr3_pssm = pssm_A, pssm_weight = pssm_weight)

  sample_chain_cdr3_multi("TRB",
                          file.path(step_dir, "TRBV_TRBJ.csv"), cdr3_baseline,
                          file.path(step_dir, "TRBV_TRBJ_cdr3.csv"),
                          n = n_cdr3_multi, len_dist = len_dist_B,
                          cdr3_pssm = pssm_B, pssm_weight = pssm_weight)

  df_paired <- pair_alpha_beta_multi(
    file.path(step_dir, "TRAV_TRAJ_cdr3.csv"),
    file.path(step_dir, "TRBV_TRBJ_cdr3.csv"),
    file.path(step_dir, "chainA_B_random_pair.csv"),
    n = n_cdr3_multi
  )
  prepare_model_csv(df_paired, peptide, mhc, file.path(step_dir, "model.csv"))

  list(
    next_model_csv = file.path(step_dir, "model.csv"),
    n_total_tcrs   = nrow(scored_tcrs),
    n_top_tcrs     = nrow(top_tcrs),
    pssm_built     = !is.null(pssm_A)
  )
}


# =============================================================================
# Main pipeline
# =============================================================================

run_motif_builder <- function(peptide, mhc,
                               base_output_dir  = BASE_OUTPUT_DIR,
                               input_dir        = INPUT_DIR,
                               cdr3_baseline,
                               known_binders_file,
                               validation_file,
                               n_iterations     = N_ITERATIONS,
                               n_enrichment     = N_ENRICHMENT,
                               n_cdr3_multi     = N_CDR3_MULTI,
                               init_perc_rank   = INIT_PERC_RANK,
                               perc_rank_mid  = PERC_RANK_MID,
                               pssm_weight      = PSSM_WEIGHT,
                               min_tcrs_pssm    = MIN_TCRS_PSSM,
                               score_col        = SCORE_COL,
                               plot_motif       = TRUE) {

  results <- list()

  # ------------------------------------------------------------------
  # Full threshold schedule: length = n_iterations
  # Step 1 always uses init_perc_rank (10), final step always 0.05,
  # middle steps follow geometric decay from perc_rank_mid.
  # e.g. N_ITERATIONS=3 → [10, perc_rank_mid, 0.05]
  # ------------------------------------------------------------------
  if (n_iterations == 1) {
    perc_rank_schedule <- init_perc_rank
  } else {
    perc_rank_schedule <- c(init_perc_rank,
                            exp(seq(log(perc_rank_mid), log(0.05),
                                    length.out = n_iterations - 1)))
  }
  message(sprintf("[%s] perc_rank schedule: %s",
                  peptide,
                  paste(sprintf("%.3f", perc_rank_schedule), collapse = " → ")))

  # ------------------------------------------------------------------
  # Build (or reload) TEMPO predictor once
  # ------------------------------------------------------------------
  predictor_rds <- file.path(BASE_OUTPUT_DIR, paste0("predictor_", peptide, ".rds"))
  predictor_rds <- build_tempo_predictor(known_binders_file, predictor_rds)

  # ------------------------------------------------------------------
  # Step 0: flat random TCR generation — no priors
  # ------------------------------------------------------------------
  message(sprintf("\n[%s] Step 0: random TCR generation (flat baseline)", peptide))

  step0_dir <- file.path(base_output_dir, sprintf("TEMPO_step0_%s", peptide))
  dir.create(step0_dir, showWarnings = FALSE, recursive = TRUE)

  generate_vj_pairs("A", input_dir, step0_dir)
  generate_vj_pairs("B", input_dir, step0_dir)
  sample_chain_cdr3("TRA", file.path(step0_dir, "TRAV_TRAJ.csv"),
                    cdr3_baseline, file.path(step0_dir, "TRAV_TRAJ_cdr3.csv"))
  sample_chain_cdr3("TRB", file.path(step0_dir, "TRBV_TRBJ.csv"),
                    cdr3_baseline, file.path(step0_dir, "TRBV_TRBJ_cdr3.csv"))
  df_paired0 <- pair_alpha_beta(
    file.path(step0_dir, "TRAV_TRAJ_cdr3.csv"),
    file.path(step0_dir, "TRBV_TRBJ_cdr3.csv"),
    file.path(step0_dir, "chainA_B_random_pair.csv")
  )
  prepare_model_csv(df_paired0, peptide, mhc, file.path(step0_dir, "model.csv"))
  current_model_csv <- file.path(step0_dir, "model.csv")

  # ------------------------------------------------------------------
  # Unified loop: N_ITERATIONS scoring+enrichment steps (steps 1..N)
  # Each step scores the current TCRs and generates the next batch.
  # ------------------------------------------------------------------
  for (iter in seq_len(n_iterations)) {
    threshold <- perc_rank_schedule[iter]
    message(sprintf("[%s] Step %d / %d (threshold = %.3f)", peptide, iter, n_iterations, threshold))

    step_result <- run_enrichment_step(
      predictor_rds     = predictor_rds,
      current_model_csv = current_model_csv,
      threshold         = threshold,
      step              = iter,
      peptide           = peptide, mhc = mhc,
      base_output_dir   = base_output_dir,
      cdr3_baseline     = cdr3_baseline,
      n_enrichment      = n_enrichment,
      n_cdr3_multi      = n_cdr3_multi,
      pssm_weight       = pssm_weight,
      min_tcrs_pssm     = min_tcrs_pssm,
      score_col         = score_col,
      plot_motif        = plot_motif
    )
    if (is.null(step_result)) break
    results[[paste0("step", iter)]] <- step_result
    current_model_csv <- step_result$next_model_csv
  }

  # ------------------------------------------------------------------
  # Final validation
  # ------------------------------------------------------------------
  message(sprintf("[%s] Final validation (AUC)", peptide))

  auc_result <- validate_motif(
    motif_file      = current_model_csv,
    validation_file = validation_file,
    output_dir      = file.path(base_output_dir, sprintf("TEMPO_validation_%s", peptide)),
    peptide         = peptide,
    mhc             = mhc
  )

  results[["final_auc"]]   <- auc_result$auc
  results[["final_auc01"]] <- auc_result$auc01
  results[["final_roc"]]   <- auc_result$roc
  results[["final_pred"]]  <- auc_result$pred

  # Flatten per-step counts for easy logging
  step_keys <- grep("^step", names(results), value = TRUE)
  results[["step_counts"]] <- do.call(rbind, lapply(step_keys, function(k) {
    data.frame(step    = as.integer(sub("^step", "", k)),
               n_total = results[[k]]$n_total_tcrs,
               n_top   = results[[k]]$n_top_tcrs)
  }))

  results
}


### ---- Run ------------------------------------------------------------------
# Guard: only execute when run directly, not when sourced by the optimizer.

baseline      <- MixTCRviz::baseline_HomoSapiens
cdr3_baseline <- baseline$countCDR3.VJL

if (!exists(".sourced_by_optimizer") || !.sourced_by_optimizer) {

  auc_summary <- list()

  for (peptide in names(dico)) {
    mhc <- dico[[peptide]]
    message(sprintf("\n========== %s / %s ==========", peptide, mhc))

    res <- run_motif_builder(
      peptide            = peptide,
      mhc                = mhc,
      cdr3_baseline      = cdr3_baseline,
      known_binders_file = file.path(peptide, paste0("A0201_", peptide, ".csv")),
      validation_file    = file.path(peptide, "validation.csv")
    )

    auc_summary[[peptide]] <- list(auc = res$final_auc, auc01 = res$final_auc01)
  }

  message("\n===== AUC Summary =====")
  for (nm in names(auc_summary)) {
    message(sprintf("  %s: AUC = %.4f  |  AUC0.1 = %.4f",
                    nm, auc_summary[[nm]]$auc, auc_summary[[nm]]$auc01))
  }

}
