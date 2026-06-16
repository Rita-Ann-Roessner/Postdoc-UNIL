library(TEMPOtrain)
library(MixTCRviz)
library(dplyr)
library(tidyr)
library(pROC)

# =============================================================================
# AF3 Motif Builder
#
# Structure-based iterative TCR motif building pipeline.
# Uses AlphaFold3 ipTM scoring with dummy partner chains.
#
# Pipeline structure:
#   Step 0:    Generate random TCRs (flat baseline) paired with dummy chain
#              → model_alpha.csv (real α + dummy β)
#              → model_beta.csv  (dummy α + real β)
#              → Run AF3 on cluster; place results as:
#                  step0/model_alpha_output.txt
#                  step0/model_beta_output.txt
#   Steps 1–N: Load AF3 scores → filter (AF3_iptm_pair_mean >= threshold)
#              → enrich V/J + CDR3 length + PSSM per chain
#              → generate new batches → run AF3 on cluster
#   "final":   Load last AF3 scores → filter → train TEMPO on top binders
#              (dummy chain columns cleared) → AUC / AUC0.1 validation
#
# Usage: set STEP before sourcing / running this script.
#   STEP = 0              → generate step-0 TCR batches
#   STEP = 1 .. N_STEPS   → enrich from step (STEP-1) AF3 scores; generate new batch
#   STEP = c(1, 2, 3)     → run multiple enrichment steps in sequence
#   STEP = "final"        → TEMPO validation using step N_STEPS AF3 scores
# =============================================================================

### ---- Configuration --------------------------------------------------------

STEP            <- c(1, 2, 3, 4)        # 0, 1, ..., N_STEPS, or "final"
N_STEPS         <- 5         # total number of enrichment steps after step 0

INPUT_DIR       <- "/Users/roessner/Documents/PostDoc/Data/MixTCRviz/data_raw/HomoSapiens"
BASE_OUTPUT_DIR <- "bck.TCR_motif_atlas"
SCORE_COL       <- "AF3_iptm_pair_mean"   # column in af3_scores_*.txt; higher = better

TOP_PCT             <- 0.01 # top % of pooled scores across all steps to keep at each step

N_PAIRS          <- 400    # V/J pairs sampled per chain from top-binder distribution
N_CDR3_MULTI     <- 3      # CDR3 sequences sampled per V/J pair (enrichment steps)
PSSM_WEIGHT      <- 0.643  # blend weight: 0 = pure baseline, 1 = pure PSSM
MIN_TCRS_PSSM    <- 30     # min top binders required to build a PSSM
LEN_DIST_COND_VJ <- FALSE  # if TRUE, sample CDR3 length | V/J; if FALSE, use marginal

PLOT_EACH_STEP     <- TRUE   # if TRUE, plot MixTCRviz motif after each enrichment step
VALIDATE_EACH_STEP <- TRUE  # if TRUE, run TEMPO validation after each enrichment step

EPITOPES_FILE <- NULL #"epitopes.txt"   # one "MHC_PEPTIDE" label per line, e.g. "A0201_ELAGIGILTV"
epitopes <- c("A0201_LLWNGPMAV", "A0201_GILGFVFTL")

if (!is.null(EPITOPES_FILE) && file.exists(EPITOPES_FILE)) {
  epitopes <- trimws(readLines(EPITOPES_FILE))
  epitopes <- epitopes[nchar(epitopes) > 0]
  message(sprintf("Loaded %d epitopes from: %s", length(epitopes), EPITOPES_FILE))
}


### ---- Dummy chain definitions ----------------------------------------------
# Most common V/J pair + polyGly CDR3, used as a neutral scaffold so AF3
# scoring reflects only the varied chain's contacts.

DUMMY_TRB <- list(
  v    = "TRBV20-1_dummy",
  j    = "TRBJ2-7",
  cdr3 = paste0("C", paste(rep("G", 12), collapse = ""), "F")
)
DUMMY_TRA <- list(
  v    = "TRAV38-2/DV8_dummy",
  j    = "TRAJ43",
  cdr3 = paste0("C", paste(rep("G", 10), collapse = ""), "F")
)


# =============================================================================
# Module 1 — Step-0 TCR generation: flat baseline sampling
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
    if (sum(p) == 0) return(NA_character_)
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
  pairs           <- expand.grid(lst[[1]], lst[[2]])
  colnames(pairs) <- genes
  write.csv(pairs, file.path(output_dir, paste0(genes[1], "_", genes[2], ".csv")),
            row.names = FALSE)
  pairs
}


sample_chain_cdr3 <- function(chain, pair_file, cdr3_baseline, output_file) {
  df           <- read.csv(pair_file)
  chain_letter <- sub("^TR", "", chain)
  v_col        <- paste0("TR", chain_letter, "V")
  j_col        <- paste0("TR", chain_letter, "J")
  df$random_CDR3 <- mapply(draw_random_cdr3,
                            chain = chain, v_seg = df[[v_col]], j_seg = df[[j_col]],
                            MoreArgs = list(cdr3_baseline = cdr3_baseline))
  df <- na.omit(df)
  write.csv(df, output_file, row.names = FALSE)
  df
}


pair_with_dummy_beta <- function(df_A, peptide, mhc_allele, species, output_file) {
  colnames(df_A)[colnames(df_A) == "random_CDR3"] <- "cdr3_TRA"
  df_A$TRBV     <- DUMMY_TRB$v
  df_A$TRBJ     <- DUMMY_TRB$j
  df_A$cdr3_TRB <- DUMMY_TRB$cdr3
  df_A$id       <- sprintf("tcr%04d", seq_len(nrow(df_A)))
  df_A$peptide  <- peptide
  df_A$MHC      <- mhc_allele
  df_A$species  <- species
  write.csv(df_A, output_file, row.names = FALSE)
  df_A
}


pair_with_dummy_alpha <- function(df_B, peptide, mhc_allele, species, output_file) {
  colnames(df_B)[colnames(df_B) == "random_CDR3"] <- "cdr3_TRB"
  df_B$TRAV     <- DUMMY_TRA$v
  df_B$TRAJ     <- DUMMY_TRA$j
  df_B$cdr3_TRA <- DUMMY_TRA$cdr3
  df_B$id       <- sprintf("tcr%04d", seq_len(nrow(df_B)))
  df_B$peptide  <- peptide
  df_B$MHC      <- mhc_allele
  df_B$species  <- species
  write.csv(df_B, output_file, row.names = FALSE)
  df_B
}


# =============================================================================
# Module 2 — PSSM from high-scoring CDR3 sequences
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
# Module 3 — Multi-CDR3 sampling with V/J + length + PSSM enrichments
# =============================================================================

draw_random_cdr3_multi <- function(chain, v_seg, j_seg, cdr3_baseline,
                                    n = 5, len_dist = NULL, len_dist_vj = NULL,
                                    cdr3_pssm = NULL, pssm_weight = 0.5) {
  key <- paste0(v_seg, "_", j_seg)
  if (!key %in% names(cdr3_baseline[[chain]])) return(character(0))

  vj_counts         <- cdr3_baseline[[chain]][[key]]
  lengths_available <- names(vj_counts)

  active_len_dist <- if (!is.null(len_dist_vj) && key %in% names(len_dist_vj)) {
    len_dist_vj[[key]]
  } else {
    len_dist
  }

  if (!is.null(active_len_dist)) {
    valid_lens <- intersect(names(active_len_dist), lengths_available)
    if (length(valid_lens) > 0) {
      probs         <- unlist(active_len_dist[valid_lens])
      probs         <- probs / sum(probs)
      selected_lens <- sample(valid_lens, size = min(n, length(valid_lens)),
                              replace = FALSE, prob = probs)
    } else {
      selected_lens <- sample(lengths_available, size = min(n, length(lengths_available)))
    }
  } else {
    selected_lens <- sample(lengths_available, size = min(n, length(lengths_available)))
  }

  seqs <- vapply(selected_lens, function(len_name) {
    count_mat <- as.matrix(vj_counts[[len_name]])
    prob_mat  <- apply(count_mat, 2, function(col) {
      if (sum(col) == 0) rep(1 / nrow(count_mat), nrow(count_mat))
      else               col / sum(col)
    })
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
                                     n = 5, len_dist = NULL, len_dist_vj = NULL,
                                     cdr3_pssm = NULL, pssm_weight = 0.5) {
  df           <- read.csv(pair_file)
  chain_letter <- sub("^TR", "", chain)
  v_col        <- paste0("TR", chain_letter, "V")
  j_col        <- paste0("TR", chain_letter, "J")

  results <- mapply(draw_random_cdr3_multi,
                    chain = chain, v_seg = df[[v_col]], j_seg = df[[j_col]],
                    MoreArgs = list(cdr3_baseline = cdr3_baseline,
                                   n           = n,
                                   len_dist    = len_dist,
                                   len_dist_vj = len_dist_vj,
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


# =============================================================================
# Module 4 — AF3 score loading and filtering
# =============================================================================

load_af3_scores <- function(scores_file, model_file, score_col = SCORE_COL) {
  scores <- read.table(scores_file, header = TRUE, sep = "\t")
  model  <- read.table(model_file,  header = TRUE, sep = "\t")

  if (!score_col %in% colnames(scores))
    stop("Column '", score_col, "' not found in AF3 scores file.\n",
         "Available: ", paste(colnames(scores), collapse = ", "))
  if (!"id" %in% colnames(scores) || !"id" %in% colnames(model))
    stop("Both the AF3 scores file and model CSV must have an 'id' column.")

  merge(model, scores[, c("id", score_col)], by = "id")
}


# =============================================================================
# Module 5 — Enrichment helpers (per-chain)
# =============================================================================

extract_vj_dist_chain <- function(top_tcrs, chain_letter) {
  v_col <- paste0("TR", chain_letter, "V")
  j_col <- paste0("TR", chain_letter, "J")
  if (!all(c(v_col, j_col) %in% colnames(top_tcrs))) return(NULL)
  keep  <- !is.na(top_tcrs[[v_col]]) & !is.na(top_tcrs[[j_col]])
  pairs <- paste(top_tcrs[[v_col]][keep], top_tcrs[[j_col]][keep], sep = "_")
  if (length(pairs) == 0) return(NULL)
  tbl <- table(pairs)
  as.list(tbl / sum(tbl))
}


sample_vj_pairs <- function(vj_dist, chain_letter, n_pairs, output_dir) {
  v_col <- paste0("TR", chain_letter, "V")
  j_col <- paste0("TR", chain_letter, "J")
  if (is.null(vj_dist))
    stop(sprintf("No V/J distribution available for chain TR%s.", chain_letter))

  sampled  <- sample(names(vj_dist), size = n_pairs, replace = TRUE,
                     prob = unlist(vj_dist))
  split_vj <- strsplit(sampled, "_")
  pairs    <- data.frame(
    V = sapply(split_vj, `[`, 1),
    J = sapply(split_vj, `[`, 2),
    stringsAsFactors = FALSE
  )
  colnames(pairs) <- c(v_col, j_col)
  write.csv(pairs, file.path(output_dir, paste0(v_col, "_", j_col, ".csv")),
            row.names = FALSE)
  pairs
}


extract_cdr3_len_dist <- function(top_tcrs, scored, chain_letter) {
  cdr3_col <- paste0("cdr3_TR", chain_letter)
  if (!cdr3_col %in% colnames(top_tcrs)) return(NULL)

  freq_top <- table(paste0("L_", nchar(top_tcrs[[cdr3_col]][!is.na(top_tcrs[[cdr3_col]])])))
  freq_all <- table(paste0("L_", nchar(scored[[cdr3_col]][!is.na(scored[[cdr3_col]])])))
  if (sum(freq_top) == 0 || sum(freq_all) == 0) return(NULL)

  freq_top <- freq_top / sum(freq_top)
  freq_all <- freq_all / sum(freq_all)
  all_lens <- union(names(freq_top), names(freq_all))
  ratio <- sapply(all_lens, function(l) {
    top_f <- if (l %in% names(freq_top)) as.numeric(freq_top[l]) else 0
    all_f <- if (l %in% names(freq_all)) as.numeric(freq_all[l]) else 0
    if (all_f == 0) return(0)
    top_f / all_f
  })
  ratio <- ratio[ratio > 0]
  if (sum(ratio) == 0) return(NULL)
  as.list(ratio / sum(ratio))
}


extract_cdr3_len_dist_vj <- function(top_tcrs, scored, chain_letter) {
  cdr3_col <- paste0("cdr3_TR", chain_letter)
  v_col    <- paste0("TR", chain_letter, "V")
  j_col    <- paste0("TR", chain_letter, "J")
  if (!all(c(cdr3_col, v_col, j_col) %in% colnames(top_tcrs))) return(NULL)

  prep <- function(df) {
    keep <- !is.na(df[[cdr3_col]]) & !is.na(df[[v_col]]) & !is.na(df[[j_col]])
    df   <- df[keep, ]
    split(paste0("L_", nchar(df[[cdr3_col]])), paste0(df[[v_col]], "_", df[[j_col]]))
  }
  top_by_vj <- prep(top_tcrs)
  all_by_vj <- prep(scored)
  if (length(top_by_vj) == 0) return(NULL)

  lapply(setNames(names(top_by_vj), names(top_by_vj)), function(vj) {
    l_top    <- top_by_vj[[vj]]
    l_all    <- if (vj %in% names(all_by_vj)) all_by_vj[[vj]] else character(0)
    freq_top <- table(l_top) / length(l_top)
    if (length(l_all) == 0) return(as.list(freq_top / sum(freq_top)))
    freq_all <- table(l_all) / length(l_all)
    all_lens <- union(names(freq_top), names(freq_all))
    ratio <- sapply(all_lens, function(l) {
      top_f <- if (l %in% names(freq_top)) as.numeric(freq_top[l]) else 0
      all_f <- if (l %in% names(freq_all)) as.numeric(freq_all[l]) else 0
      if (all_f == 0) return(0)
      top_f / all_f
    })
    ratio <- ratio[ratio > 0]
    if (sum(ratio) == 0) return(NULL)
    as.list(ratio / sum(ratio))
  })
}


# Pool all raw scored TCRs from steps 0 through max_step (inclusive).
load_all_scored <- function(chain_letter, max_step, label, base_output_dir) {
  chain_name <- if (chain_letter == "A") "alpha" else "beta"
  dfs <- lapply(0:max_step, function(k) {
    scores_file <- file.path(base_output_dir, label, sprintf("step%d", k),
                             sprintf("model_%s_output.txt", chain_name))
    model_file  <- file.path(base_output_dir, label, sprintf("step%d", k),
                             sprintf("model_%s_input.txt",  chain_name))
    if (!file.exists(scores_file) || !file.exists(model_file)) return(NULL)
    load_af3_scores(scores_file, model_file)
  })
  dfs <- Filter(Negate(is.null), dfs)
  if (length(dfs) == 0) return(NULL)
  bind_rows(dfs)
}


# =============================================================================
# Module 6 — TEMPO final validation
# Trains TEMPO directly from the top-binder motif file (no percentile-rank
# predictor needed); uses raw score for ROC.
# =============================================================================

validate_motif <- function(motif_file, validation_file, output_dir,
                            peptide, mhc, score_col = "score", plot = TRUE) {
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
  if (length(pred_files) == 0) stop("No validation output found in ", output_dir)
  pred_df <- read.csv(pred_files[1])

  if (!"Label" %in% colnames(pred_df))
    stop("'Label' column not found in validation output; check that validation.csv has 'Label'.")
  if (!score_col %in% colnames(pred_df))
    stop("Column '", score_col, "' not found in validation output.\n",
         "Available: ", paste(colnames(pred_df), collapse = ", "))

  pred_df   <- pred_df[!is.na(pred_df[[score_col]]), ]
  roc_obj   <- roc(pred_df$Label, pred_df[[score_col]], quiet = TRUE)
  auc_val   <- as.numeric(auc(roc_obj))
  auc01_val <- as.numeric(auc(roc_obj,
                               partial.auc         = c(1, 0.9),
                               partial.auc.correct = TRUE,
                               partial.auc.focus   = "specificity"))

  message(sprintf("  AUC = %.4f  |  AUC0.1 = %.4f", auc_val, auc01_val))

  if (plot) {
    plot(roc_obj,
         main = sprintf("ROC - %s_%s", mhc, peptide),
         col  = "#2166ac", lwd = 2, xaxt = "n",
         xlab = "FPR (1 - Specificity)", ylab = "TPR (Sensitivity)")
    axis(1, at = seq(1, 0, -0.2), labels = seq(0, 1, 0.2))
    legend("bottomright", bty = "o", bg = "white", box.col = "grey80",
           legend = c(sprintf("AUC    = %.4f", auc_val),
                      sprintf("AUC0.1 = %.4f", auc01_val)))
  }

  list(auc = auc_val, auc01 = auc01_val, roc = roc_obj, pred = pred_df)
}


# =============================================================================
# Module 7 — Per-step motif plot and validation helpers
# =============================================================================

# Plot a MixTCRviz motif for one chain's top binders vs the default baseline.
plot_step_motif <- function(top_tcrs, chain_letter, label, step, step_dir) {
  chain_name <- if (chain_letter == "A") "alpha" else "beta"
  out_dir    <- file.path(step_dir, sprintf("motif_%s", chain_name))
  n_top      <- nrow(top_tcrs)
  message(sprintf("  Plotting motif (step %d, chain %s, n_top=%d)", step, chain_letter, n_top))
  MixTCRviz(
    input1      = top_tcrs,
    output.path = out_dir,
    plot        = TRUE,
    renormVJ    = TRUE,
    set.title   = sprintf("%s - step %d %s top binders (n=%d)", label, step, chain_name, n_top),
    verbose     = 0
  )
}


# Combine alpha and beta top binders into a single TEMPO training df,
# clearing the dummy-chain columns in each.
build_motif_df <- function(top_alpha, top_beta, model_label) {
  motif_rows <- list()
  if (!is.null(top_alpha) && nrow(top_alpha) > 0) {
    df_a <- top_alpha[, c("TRAV", "TRAJ", "cdr3_TRA", "TRBV", "TRBJ", "cdr3_TRB"), drop = FALSE]
    df_a$TRBV <- df_a$TRBJ <- df_a$cdr3_TRB <- NA_character_
    df_a$model <- model_label
    motif_rows[["alpha"]] <- df_a
  }
  if (!is.null(top_beta) && nrow(top_beta) > 0) {
    df_b <- top_beta[, c("TRAV", "TRAJ", "cdr3_TRA", "TRBV", "TRBJ", "cdr3_TRB"), drop = FALSE]
    df_b$TRAV <- df_b$TRAJ <- df_b$cdr3_TRA <- NA_character_
    df_b$model <- model_label
    motif_rows[["beta"]] <- df_b
  }
  if (length(motif_rows) == 0) return(NULL)
  do.call(rbind, motif_rows)
}


# =============================================================================
# Step-level pipeline functions
# =============================================================================

run_step0 <- function(peptide, mhc_allele, label, cdr3_baseline,
                       base_output_dir, input_dir, species = "HomoSapiens") {
  step_dir <- file.path(base_output_dir, label, "step0")
  dir.create(step_dir, showWarnings = FALSE, recursive = TRUE)

  message(sprintf("[%s] Step 0: generating flat-random TCR batches with dummy partner chains", label))

  # Alpha chain: real TRA + dummy TRB
  generate_vj_pairs("A", input_dir, step_dir)
  df_A <- sample_chain_cdr3("TRA",
                              file.path(step_dir, "TRAV_TRAJ.csv"),
                              cdr3_baseline,
                              file.path(step_dir, "TRAV_TRAJ_cdr3.csv"))
  pair_with_dummy_beta(df_A, peptide, mhc_allele, species,
                        file.path(step_dir, "model_alpha.csv"))

  # Beta chain: dummy TRA + real TRB
  generate_vj_pairs("B", input_dir, step_dir)
  df_B <- sample_chain_cdr3("TRB",
                              file.path(step_dir, "TRBV_TRBJ.csv"),
                              cdr3_baseline,
                              file.path(step_dir, "TRBV_TRBJ_cdr3.csv"))
  pair_with_dummy_alpha(df_B, peptide, mhc_allele, species,
                         file.path(step_dir, "model_beta.csv"))

  n_alpha <- nrow(read.csv(file.path(step_dir, "model_alpha.csv")))
  n_beta  <- nrow(read.csv(file.path(step_dir, "model_beta.csv")))
  message(sprintf("[%s] Step 0 done: %d alpha-batch TCRs, %d beta-batch TCRs.", label, n_alpha, n_beta))
  message(sprintf("[%s] → Run AF3 on cluster, then place results as:\n  %s\n  %s",
                  label,
                  file.path(step_dir, "model_alpha_output.txt"),
                  file.path(step_dir, "model_beta_output.txt")))

  invisible(step_dir)
}


# Enrich one chain from previous AF3 scores and generate a new model CSV.
# Returns top_tcrs (filtered) or NULL if too few binders.
enrich_one_chain <- function(chain_letter, step, label, peptide, mhc_allele, species,
                              base_output_dir, cdr3_baseline,
                              n_pairs, n_cdr3_multi, pssm_weight, min_tcrs_pssm,
                              top_pct, len_dist_conditioned_on_vj = LEN_DIST_COND_VJ) {

  chain_name    <- if (chain_letter == "A") "alpha" else "beta"
  prev_step_dir <- file.path(base_output_dir, label, sprintf("step%d", step - 1))
  step_dir      <- file.path(base_output_dir, label, sprintf("step%d", step))
  dir.create(step_dir, showWarnings = FALSE, recursive = TRUE)

  scores_file <- file.path(prev_step_dir, sprintf("model_%s_output.txt", chain_name))
  model_file  <- file.path(prev_step_dir, sprintf("model_%s_input.txt",  chain_name))

  if (!file.exists(scores_file))
    stop(sprintf("[%s] AF3 scores file missing: %s\nRun AF3 on the cluster first.", label, scores_file))
  if (!file.exists(model_file))
    stop(sprintf("[%s] Model input file missing: %s\nThis file is produced by Script_build_struct_inputs_mod.R on the cluster.", label, model_file))

  # Pool all scored TCRs from steps 0 through step-1, then take top X%
  all_scored <- load_all_scored(chain_letter, step - 1, label, base_output_dir)
  n_keep     <- max(1L, ceiling(nrow(all_scored) * top_pct))
  top_tcrs   <- head(all_scored[order(-all_scored[[SCORE_COL]]), ], n_keep)
  message(sprintf("  [chain %s] top %.0f%% = %d / %d pooled TCRs (min score = %.3f)",
                  chain_letter, top_pct * 100, nrow(top_tcrs), nrow(all_scored),
                  min(top_tcrs[[SCORE_COL]])))

  write.csv(top_tcrs, file.path(prev_step_dir, sprintf("top_tcrs_%s.csv", chain_name)),
            row.names = FALSE)

  if (nrow(top_tcrs) < 10) {
    warning(sprintf("[%s] Too few top binders for chain %s (n=%d) — skipping enrichment.",
                    label, chain_letter, nrow(top_tcrs)))
    return(NULL)
  }

  # (1) V/J distribution from top binders
  vj_dist <- extract_vj_dist_chain(top_tcrs, chain_letter)
  sample_vj_pairs(vj_dist, chain_letter, n_pairs, step_dir)

  # (2) CDR3 length distribution (relative enrichment: top binders vs all scored)
  if (len_dist_conditioned_on_vj) {
    len_dist    <- NULL
    len_dist_vj <- extract_cdr3_len_dist_vj(top_tcrs, all_scored, chain_letter)
  } else {
    len_dist    <- extract_cdr3_len_dist(top_tcrs, all_scored, chain_letter)
    len_dist_vj <- NULL
  }

  # (3) PSSM from top binders' CDR3s
  cdr3_col <- paste0("cdr3_TR", chain_letter)
  if (nrow(top_tcrs) >= min_tcrs_pssm) {
    cdr3_pssm <- build_cdr3_pssm(top_tcrs[[cdr3_col]])
    message(sprintf("  PSSM built from %d top binders (chain %s)", nrow(top_tcrs), chain_letter))
  } else {
    cdr3_pssm <- NULL
    message(sprintf("  Too few for PSSM (n=%d < %d, chain %s)", nrow(top_tcrs), min_tcrs_pssm, chain_letter))
  }

  # Generate enriched CDR3 batch
  chain_tr  <- paste0("TR", chain_letter)
  pair_file <- file.path(step_dir, sprintf("TR%sV_TR%sJ.csv", chain_letter, chain_letter))
  cdr3_file <- file.path(step_dir, sprintf("TR%sV_TR%sJ_cdr3.csv", chain_letter, chain_letter))
  sample_chain_cdr3_multi(chain_tr, pair_file, cdr3_baseline, cdr3_file,
                           n = n_cdr3_multi, len_dist = len_dist, len_dist_vj = len_dist_vj,
                           cdr3_pssm = cdr3_pssm, pssm_weight = pssm_weight)

  df_chain  <- read.csv(cdr3_file)
  model_out <- file.path(step_dir, sprintf("model_%s.csv", chain_name))
  if (chain_letter == "A") {
    pair_with_dummy_beta(df_chain, peptide, mhc_allele, species, model_out)
  } else {
    pair_with_dummy_alpha(df_chain, peptide, mhc_allele, species, model_out)
  }

  top_tcrs
}


run_enrich_step <- function(step, peptide, mhc_allele, label,
                             cdr3_baseline, base_output_dir,
                             n_pairs, n_cdr3_multi, pssm_weight, min_tcrs_pssm,
                             top_pct, species = "HomoSapiens",
                             len_dist_conditioned_on_vj = LEN_DIST_COND_VJ,
                             plot_each_step     = PLOT_EACH_STEP,
                             validate_each_step = VALIDATE_EACH_STEP,
                             validation_file    = NULL, mhc = NULL) {

  prev_step_dir <- file.path(base_output_dir, label, sprintf("step%d", step - 1))
  step_dir      <- file.path(base_output_dir, label, sprintf("step%d", step))

  message(sprintf("[%s] Step %d: enrichment + generation (top %.0f%% pooled across all steps)",
                  label, step, top_pct * 100))

  top_alpha <- enrich_one_chain("A", step, label, peptide, mhc_allele, species,
                                base_output_dir, cdr3_baseline,
                                n_pairs, n_cdr3_multi, pssm_weight, min_tcrs_pssm,
                                top_pct, len_dist_conditioned_on_vj)
  top_beta  <- enrich_one_chain("B", step, label, peptide, mhc_allele, species,
                                base_output_dir, cdr3_baseline,
                                n_pairs, n_cdr3_multi, pssm_weight, min_tcrs_pssm,
                                top_pct, len_dist_conditioned_on_vj)

  # top_alpha/top_beta are already the pooled top X% — use directly for motif and training
  if (plot_each_step) {
    if (!is.null(top_alpha) && nrow(top_alpha) >= 10)
      plot_step_motif(top_alpha, "A", label, step - 1, prev_step_dir)
    if (!is.null(top_beta) && nrow(top_beta) >= 10)
      plot_step_motif(top_beta, "B", label, step - 1, prev_step_dir)
  }

  model_label <- paste0(mhc, "_", peptide)
  motif_df    <- build_motif_df(top_alpha, top_beta, model_label)
  motif_file  <- NULL
  if (!is.null(motif_df)) {
    motif_file <- file.path(prev_step_dir, "motif_top_binders.csv")
    write.csv(motif_df, motif_file, row.names = FALSE)
  }

  if (validate_each_step && !is.null(motif_file) && !is.null(validation_file) && !is.null(mhc)) {
    validate_motif(
      motif_file      = motif_file,
      validation_file = validation_file,
      output_dir      = file.path(prev_step_dir, "validation"),
      peptide         = peptide,
      mhc             = mhc
    )
  }

  message(sprintf("[%s] Step %d done.", label, step))
  message(sprintf("[%s] → Run AF3 on cluster, then place results as:\n  %s\n  %s",
                  label,
                  file.path(step_dir, "model_alpha_output.txt"),
                  file.path(step_dir, "model_beta_output.txt")))

  invisible(list(top_alpha = top_alpha, top_beta = top_beta))
}


run_final_validation <- function(label, peptide, mhc, mhc_allele,
                                  base_output_dir, n_steps, top_pct,
                                  validation_file) {
  message(sprintf("[%s] Final: pooling scores from all %d steps (top %.0f%%)",
                  label, n_steps, top_pct * 100))

  collect_chain <- function(chain_letter) {
    all_scored <- load_all_scored(chain_letter, n_steps, label, base_output_dir)
    if (is.null(all_scored)) {
      warning(sprintf("[%s] No scored TCRs found for chain %s.", label, chain_letter))
      return(NULL)
    }
    n_keep   <- max(1L, ceiling(nrow(all_scored) * top_pct))
    top_tcrs <- head(all_scored[order(-all_scored[[SCORE_COL]]), ], n_keep)
    message(sprintf("  [chain %s] top %.0f%% = %d / %d pooled TCRs (min score = %.3f)",
                    chain_letter, top_pct * 100, nrow(top_tcrs), nrow(all_scored),
                    min(top_tcrs[[SCORE_COL]])))
    top_tcrs
  }

  top_alpha <- collect_chain("A")
  top_beta  <- collect_chain("B")

  model_label <- paste0(mhc, "_", peptide)
  motif_df    <- build_motif_df(top_alpha, top_beta, model_label)

  if (is.null(motif_df)) {
    warning(sprintf("[%s] No top binders found — cannot run TEMPO validation.", label))
    return(NULL)
  }
  motif_file <- file.path(base_output_dir, label, "motif_top_binders.csv")
  write.csv(motif_df, motif_file, row.names = FALSE)
  message(sprintf("[%s] Motif training set: %d TCRs written to %s",
                  label, nrow(motif_df), motif_file))

  val_dir    <- file.path(base_output_dir, label, "validation")
  auc_result <- validate_motif(
    motif_file      = motif_file,
    validation_file = validation_file,
    output_dir      = val_dir,
    peptide         = peptide,
    mhc             = mhc
  )

  list(auc = auc_result$auc, auc01 = auc_result$auc01)
}


# =============================================================================
# Main
# =============================================================================

baseline      <- MixTCRviz::baseline_HomoSapiens
cdr3_baseline <- baseline$countCDR3.VJL

auc_summary <- list()

for (epitope in epitopes) {
  mhc        <- sub("_.*", "", epitope)                                         # e.g. "A0201"
  peptide    <- sub("^[^_]*_", "", epitope)                                     # e.g. "ELAGIGILTV"
  label      <- epitope                                                          # e.g. "A0201_ELAGIGILTV"
  mhc_allele <- if (grepl("^HLA_", mhc)) mhc else paste0("HLA_", mhc)         # e.g. "HLA_A0201"

  message(sprintf("\n========== %s (STEP = %s) ==========", epitope, STEP))

  if (identical(STEP, 0) || identical(STEP, 0L)) {

    run_step0(
      peptide         = peptide,
      mhc_allele      = mhc_allele,
      label           = label,
      cdr3_baseline   = cdr3_baseline,
      base_output_dir = BASE_OUTPUT_DIR,
      input_dir       = INPUT_DIR
    )

  } else if (is.numeric(STEP) && all(STEP >= 1 & STEP <= N_STEPS)) {

    for (s in sort(unique(STEP))) {
      run_enrich_step(
        step                       = s,
        peptide                    = peptide,
        mhc_allele                 = mhc_allele,
        label                      = label,
        cdr3_baseline              = cdr3_baseline,
        base_output_dir            = BASE_OUTPUT_DIR,
        n_pairs                    = N_PAIRS,
        n_cdr3_multi               = N_CDR3_MULTI,
        pssm_weight                = PSSM_WEIGHT,
        min_tcrs_pssm              = MIN_TCRS_PSSM,
        top_pct                    = TOP_PCT,
        len_dist_conditioned_on_vj = LEN_DIST_COND_VJ,
        plot_each_step             = PLOT_EACH_STEP,
        validate_each_step         = VALIDATE_EACH_STEP,
        validation_file            = file.path(BASE_OUTPUT_DIR, label, "validation.csv"),
        mhc                        = mhc
      )
    }

  } else if (identical(STEP, "final")) {

    res <- run_final_validation(
      label           = label,
      peptide         = peptide,
      mhc             = mhc,
      mhc_allele      = mhc_allele,
      base_output_dir = BASE_OUTPUT_DIR,
      n_steps         = N_STEPS,
      top_pct         = TOP_PCT,
      validation_file = file.path(BASE_OUTPUT_DIR, label, "validation.csv")
    )
    if (!is.null(res)) auc_summary[[epitope]] <- res

  } else {
    stop(sprintf("Invalid STEP = %s. Must be 0, 1..%d, or \"final\".", STEP, N_STEPS))
  }
}

if (identical(STEP, "final") && length(auc_summary) > 0) {
  auc_df <- data.frame(
    epitope = names(auc_summary),
    auc     = sapply(auc_summary, `[[`, "auc"),
    auc01   = sapply(auc_summary, `[[`, "auc01"),
    row.names = NULL
  )
  write.csv(auc_df, file.path(BASE_OUTPUT_DIR, "auc_summary.csv"), row.names = FALSE)
  message("\n===== AUC Summary =====")
  for (nm in names(auc_summary)) {
    message(sprintf("  %s: AUC = %.4f  |  AUC0.1 = %.4f",
                    nm, auc_summary[[nm]]$auc, auc_summary[[nm]]$auc01))
  }
  message(sprintf("  Saved to: %s", file.path(BASE_OUTPUT_DIR, "auc_summary.csv")))
}
