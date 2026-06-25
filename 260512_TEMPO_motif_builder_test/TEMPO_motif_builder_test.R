library(TEMPOtrain)
library(MixTCRviz)
library(dplyr)
library(tidyr)
library(pROC)

# =============================================================================
# TEMPO Motif Builder
#
# Pipeline structure:
#   Step 0:    Generate random TCRs (flat baseline) → score (INIT_PERC_RANK)
#              → plot motif → top_tcrs
#   Steps 1–N: Enrich from top_tcrs (V/J, CDR3 length, PSSM) → generate
#              new TCRs → score (PERC_RANK) → plot motif → top_tcrs
#   Final:     Validate last model on labelled validation.csv → AUC
#
# NOTE: set SCORE_COL to the column name TEMPOtrain writes in its prediction
#       output (check head() of the pred file after the first run).
# =============================================================================

### ---- Configuration --------------------------------------------------------

INPUT_DIR       <- "/Users/roessner/Documents/PostDoc/Data/MixTCRviz/data_raw/HomoSapiens"
BASE_OUTPUT_DIR <- "test"
SCORE_COL       <- "perc_rank" # lower = better binder (HLA convention)
N_ITERATIONS    <- 3           # scoring+enrichment steps (steps 1–N); step 0 is always flat random
N_PAIRS         <- 400         # V/J pairs sampled per chain from top-binder distribution
N_CDR3_MULTI    <- 3           # CDR3 sequences sampled per V/J pair (iterations > 0)
INIT_PERC_RANK  <- 10          # threshold for step 0; decays each step by DECAY_FACTOR
DECAY_FACTOR    <- 0.135         # multiplicative decay per step: step k uses INIT × DECAY^k
PSSM_WEIGHT     <- 1.0        # blend weight: 0 = pure baseline, 1 = pure PSSM
MIN_TCRS_PSSM   <- 30         # minimum high-scorers required before using a PSSM
LEN_DIST_COND_VJ <- FALSE      # if TRUE, sample CDR3 length | V/J pair; if FALSE, use marginal length dist
VJ_PRIOR_STRENGTH <- 20        # alpha: repertoire-prior pseudocount weight for V/J shrinkage
                                # (in evidence-count units; higher = more shrinkage toward baseline)

EPITOPES_FILE <- NULL #"IMMREP23/epitopes.txt"   # path to a plain-text file with one epitope per line,
                        # e.g. "epitopes.txt"; set to NULL to use the list below
epitopes <- c(
  "A0201_LLWNGPMAV",
  "A0201_ELAGIGILTV",
  "A0201_GILGFVFTL"
)

if (!is.null(EPITOPES_FILE) && file.exists(EPITOPES_FILE)) {
  epitopes <- trimws(readLines(EPITOPES_FILE))
  epitopes <- epitopes[nchar(epitopes) > 0]   # drop blank lines
  message(sprintf("Loaded %d epitopes from: %s", length(epitopes), EPITOPES_FILE))
}


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
                                    n              = 5,
                                    len_dist       = NULL,
                                    len_dist_vj    = NULL,
                                    cdr3_pssm      = NULL,
                                    pssm_weight    = 0.5) {
  key <- paste0(v_seg, "_", j_seg)
  if (!key %in% names(cdr3_baseline[[chain]])) return(character(0))

  vj_counts         <- cdr3_baseline[[chain]][[key]]
  lengths_available <- names(vj_counts)

  # (2) Sample n lengths with replacement from enriched distribution.
  # Priority: VJ-conditioned > marginal > uniform from baseline.
  active_len_dist <- if (!is.null(len_dist_vj) && key %in% names(len_dist_vj)) {
    len_dist_vj[[key]]   # P(Length | V, J)
  } else {
    len_dist             # P(Length) marginal — or NULL
  }

  if (!is.null(active_len_dist)) {
    valid_lens <- intersect(names(active_len_dist), lengths_available)
    if (length(valid_lens) > 0) {
      probs         <- unlist(active_len_dist[valid_lens])
      probs         <- probs / sum(probs)
      selected_lens <- sample(valid_lens, size = min(n, length(valid_lens)),
                              replace = FALSE, prob = probs)
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

    # (3) Blend with PSSM when available and dimensions match.
    # Per-position weight (IC-modulated): trust the PSSM more where the
    # baseline is uninformative (low IC — the variable junction) and less
    # where the baseline is conserved (high IC — the anchors). IC is computed
    # from the normalized baseline column, so the weight is scale-free
    # (independent of baseline depth M and of the number of top binders).
    #   w(pos) = pssm_weight * (1 - IC_baseline(pos) / log2(K))
    # with IC_baseline(pos) = log2(K) - H(pos), K = alphabet size.
    if (!is.null(cdr3_pssm) && !is.null(cdr3_pssm[[len_name]])) {
      pssm      <- cdr3_pssm[[len_name]]
      aa_common <- intersect(rownames(prob_mat), rownames(pssm))
      if (length(aa_common) >= 10 && ncol(pssm) == ncol(prob_mat)) {
        max_ic <- log2(nrow(prob_mat))                 # log2(K); K = 20 for AAs
        ic_pos <- apply(prob_mat, 2, function(col) {
          p <- col[col > 0]
          max_ic + sum(p * log2(p))                    # = max_ic - H(col)
        })
        w_pos <- pssm_weight * (1 - ic_pos / max_ic)   # per-position PSSM weight
        w_pos <- pmax(0, pmin(pssm_weight, w_pos))     # numerical safety
        for (pos in seq_len(ncol(prob_mat))) {
          prob_mat[aa_common, pos] <- (1 - w_pos[pos]) * prob_mat[aa_common, pos] +
                                            w_pos[pos]  * pssm[aa_common, pos]
        }
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
                                     len_dist_vj = NULL,
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


pair_alpha_beta_multi <- function(file_a, file_b, output_file) {
  df_A <- read.csv(file_a)
  df_B <- read.csv(file_b)
  colnames(df_A)[colnames(df_A) == "random_CDR3"] <- "cdr3_TRA"
  colnames(df_B)[colnames(df_B) == "random_CDR3"] <- "cdr3_TRB"

  n_A <- nrow(df_A); n_B <- nrow(df_B); n_max <- max(n_A, n_B)

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

plot_iter_motif <- function(top_tcrs, label, iter, output_dir, all_tcrs = NULL) {
  baseline_label <- if (is.null(all_tcrs)) "standard baseline" else "iter TCRs"
  message(sprintf("  Plotting MixTCRviz motif (iter %d, n_top=%d, baseline=%s)",
                  iter, nrow(top_tcrs), baseline_label))
  args <- list(
    input1      = top_tcrs,
    output.path = output_dir,
    plot        = TRUE,
    renormVJ    = TRUE,
    set.title   = sprintf("%s - iter %d top binders (n=%d)", label, iter, nrow(top_tcrs)),
    verbose     = 0
  )
  if (!is.null(all_tcrs)) args$input2 <- all_tcrs
  do.call(MixTCRviz, args)
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

# Repertoire V/J prior P_baseline(V, J), straight from MixTCRviz's
# precomputed joint V/J frequency matrix (already normalized to sum to 1).
# Returns a named list over "TRAV_TRAJ" pair strings.
extract_vj_baseline_prior <- function(chain) {
  chain_full <- paste0("TR", chain)
  m <- baseline_HomoSapiens$countVJ[[chain_full]]
  v <- as.vector(m)
  names(v) <- outer(rownames(m), colnames(m), paste, sep = "_")
  as.list(v[v > 0])
}



# Evidence-weighted, enrichment-based, baseline-shrunk V/J distribution.
#   evidence(VJ)   = max(0, n_top(VJ) - n_all(VJ) * p_global)   [excess over
#                    background pass rate — a VJ passing at exactly the
#                    background rate contributes ~0 evidence, so it gets no
#                    boost from raw passer counts alone]
#   posterior(VJ) ∝ alpha * P_baseline(VJ) + evidence(VJ)        [Dirichlet/
#                    empirical-Bayes shrinkage toward the repertoire prior,
#                    weighted by how much evidence has accumulated]
# Restricted to V/J pairs present in vj_baseline_prior, since those are the
# only ones we can later draw CDR3 sequences for.
extract_vj_dist_shrunk <- function(top_tcrs, scored, chain, vj_baseline_prior, alpha) {
  v_col <- paste0("TR", chain, "V")
  j_col <- paste0("TR", chain, "J")
  if (!all(c(v_col, j_col) %in% colnames(top_tcrs))) return(NULL)

  key_of <- function(df) {
    keep <- !is.na(df[[v_col]]) & !is.na(df[[j_col]])
    paste(df[[v_col]][keep], df[[j_col]][keep], sep = "_")
  }
  top_keys <- key_of(top_tcrs)
  all_keys <- key_of(scored)
  if (length(all_keys) == 0) return(NULL)

  n_top    <- table(top_keys)
  n_all    <- table(all_keys)
  p_global <- length(top_keys) / length(all_keys)

  vj_pairs <- names(vj_baseline_prior)
  evidence <- sapply(vj_pairs, function(vj) {
    nt <- if (vj %in% names(n_top)) n_top[[vj]] else 0
    na <- if (vj %in% names(n_all)) n_all[[vj]] else 0
    max(0, nt - na * p_global)
  })

  prior     <- unlist(vj_baseline_prior[vj_pairs])
  posterior <- alpha * prior + evidence
  if (sum(posterior) == 0) return(NULL)
  as.list(posterior / sum(posterior))
}


# Sample n_pairs (V, J) pairs per chain from the joint top-binder distribution.
# Saves a pair CSV in output_dir (same format as generate_vj_pairs).
sample_vj_pairs <- function(vj_dist, chain, n_pairs, output_dir) {
  v_col <- paste0("TR", chain, "V")
  j_col <- paste0("TR", chain, "J")

  dist <- vj_dist[[chain]]
  if (is.null(dist))
    stop(sprintf("No joint V/J distribution for chain TR%s", chain))

  sampled    <- sample(names(dist), size = n_pairs, replace = TRUE,
                       prob = unlist(dist))
  split_vj   <- strsplit(sampled, "_")
  v_genes    <- sapply(split_vj, `[`, 1)
  j_genes    <- sapply(split_vj, `[`, 2)

  pairs <- data.frame(v_genes, j_genes, stringsAsFactors = FALSE)
  colnames(pairs) <- c(v_col, j_col)
  write.csv(pairs, file.path(output_dir, paste0(v_col, "_", j_col, ".csv")),
            row.names = FALSE)
  pairs
}


extract_cdr3_len_dist <- function(top_tcrs, scored, chain) {
  cdr3_col <- paste0("cdr3_TR", chain)
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


# Per-(V,J) CDR3 length distribution from top binders (conditioned on VJ).
# Returns a named list keyed by "TRAV_TRAJ" (or "TRBV_TRBJ") strings.
# Each element is a length-probability list (same format as extract_cdr3_len_dist).
# Falls back gracefully: callers should use the marginal len_dist when the VJ
# key is absent (e.g. rare pairs with zero top-binder observations).
extract_cdr3_len_dist_vj <- function(top_tcrs, scored, chain) {
  cdr3_col <- paste0("cdr3_TR", chain)
  v_col    <- paste0("TR", chain, "V")
  j_col    <- paste0("TR", chain, "J")
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



# =============================================================================
# Module 6 — AUC validation
# Train TEMPO on the final motif; score validation.csv; compute AUC.
# =============================================================================

validate_motif <- function(motif_file, validation_file, output_dir,
                            peptide, mhc,
                            score_col = "score", plot = TRUE) {
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
         col  = "#2166ac", lwd = 2,
         xaxt = "n",
         xlab = "FPR (1 - Specificity)",
         ylab = "TPR (Sensitivity)")
    axis(1, at = seq(1, 0, -0.2), labels = seq(0, 1, 0.2))
    legend("bottomright", bty = "o", bg = "white",
           box.col = "grey80",
           legend  = c(sprintf("AUC    = %.4f", auc_val),
                       sprintf("AUC0.1 = %.4f", auc01_val)))
  }

  list(auc = auc_val, auc01 = auc01_val, roc = roc_obj, pred = pred_df)
}


# =============================================================================
# Module 7a — Score a model.csv and filter top binders
#
# Scores model_csv with the cached TEMPO predictor, filters rows whose
# score_col < threshold, optionally plots a MixTCRviz motif.
# Returns a list with top_tcrs, all_tcrs (scored), and counts.
# =============================================================================

score_step <- function(predictor_rds, model_csv, threshold,
                       step, label,
                       base_output_dir, score_col,
                       plot_motif = TRUE) {

  step_dir <- file.path(base_output_dir, label, sprintf("step%d", step))
  dir.create(step_dir, showWarnings = FALSE, recursive = TRUE)

  pred_output <- run_tempo_scoring(
    predictor_rds = predictor_rds,
    pred_file     = model_csv,
    output_dir    = file.path(step_dir, "TEMPO_scoring")
  )
  scored_tcrs <- load_tempo_scores(pred_output, model_csv, score_col)
  top_tcrs    <- scored_tcrs[scored_tcrs[[score_col]] < threshold, ]

  message(sprintf("  [Step %d] %d / %d TCRs with %s < %.2f",
                  step, nrow(top_tcrs), nrow(scored_tcrs), score_col, threshold))

  write.csv(top_tcrs, file.path(step_dir, "top_tcrs.csv"), row.names = FALSE)

  if (plot_motif && nrow(top_tcrs) >= 10) {
    plot_iter_motif(top_tcrs, label, step,
                    output_dir = file.path(step_dir, "motif"),
                    all_tcrs   = scored_tcrs)
  }

  list(
    top_tcrs = top_tcrs,
    all_tcrs = scored_tcrs,
    n_total  = nrow(scored_tcrs),
    n_top    = nrow(top_tcrs)
  )
}


# =============================================================================
# Module 7b — Enrich from top binders and generate the next TCR batch
#
# Applies V/J enrichment, CDR3 length re-weighting, and PSSM blending
# derived from top_tcrs / all_tcrs, then generates a new model.csv.
# Returns the path to the new model.csv, or NULL if too few top binders.
# =============================================================================

enrich_generate_step <- function(top_tcrs, all_tcrs,
                                  step, label, peptide, mhc,
                                  base_output_dir, cdr3_baseline,
                                  n_pairs, n_cdr3_multi,
                                  pssm_weight, min_tcrs_pssm,
                                  vj_baseline_prior, vj_prior_strength,
                                  len_dist_conditioned_on_vj = LEN_DIST_COND_VJ) {

  if (nrow(top_tcrs) < 10) {
    warning(sprintf("Too few high-scoring TCRs for step %d enrichment — stopping.", step))
    return(NULL)
  }

  step_dir <- file.path(base_output_dir, label, sprintf("step%d", step))
  dir.create(step_dir, showWarnings = FALSE, recursive = TRUE)

  # (1) V/J distribution: enrichment evidence shrunk toward repertoire baseline
  vj_dist <- list(
    A = extract_vj_dist_shrunk(top_tcrs, all_tcrs, "A", vj_baseline_prior$A, vj_prior_strength),
    B = extract_vj_dist_shrunk(top_tcrs, all_tcrs, "B", vj_baseline_prior$B, vj_prior_strength)
  )
  sample_vj_pairs(vj_dist, "A", n_pairs, step_dir)
  sample_vj_pairs(vj_dist, "B", n_pairs, step_dir)

  # (2) CDR3 length distribution from top binders
  # Either conditioned on V/J pair, or marginal across all top binders.
  if (len_dist_conditioned_on_vj) {
    len_dist_A    <- NULL   # not used when conditioning on VJ
    len_dist_B    <- NULL
    len_dist_vj_A <- extract_cdr3_len_dist_vj(top_tcrs, all_tcrs, "A")
    len_dist_vj_B <- extract_cdr3_len_dist_vj(top_tcrs, all_tcrs, "B")
  } else {
    len_dist_A    <- extract_cdr3_len_dist(top_tcrs, all_tcrs, "A")
    len_dist_B    <- extract_cdr3_len_dist(top_tcrs, all_tcrs, "B")
    len_dist_vj_A <- NULL
    len_dist_vj_B <- NULL
  }

  # (3) PSSM from top binders
  if (nrow(top_tcrs) >= min_tcrs_pssm) {
    pssm_A <- build_cdr3_pssm(top_tcrs$cdr3_TRA)
    pssm_B <- build_cdr3_pssm(top_tcrs$cdr3_TRB)
    message(sprintf("  PSSM built from %d high-scoring TCRs", nrow(top_tcrs)))
  } else {
    pssm_A <- pssm_B <- NULL
    message(sprintf("  Too few high-scorers for PSSM (%d < %d)", nrow(top_tcrs), min_tcrs_pssm))
  }

  # Generate enriched TCR batch
  sample_chain_cdr3_multi("TRA",
                          file.path(step_dir, "TRAV_TRAJ.csv"), cdr3_baseline,
                          file.path(step_dir, "TRAV_TRAJ_cdr3.csv"),
                          n = n_cdr3_multi, len_dist = len_dist_A,
                          len_dist_vj = len_dist_vj_A,
                          cdr3_pssm = pssm_A, pssm_weight = pssm_weight)

  sample_chain_cdr3_multi("TRB",
                          file.path(step_dir, "TRBV_TRBJ.csv"), cdr3_baseline,
                          file.path(step_dir, "TRBV_TRBJ_cdr3.csv"),
                          n = n_cdr3_multi, len_dist = len_dist_B,
                          len_dist_vj = len_dist_vj_B,
                          cdr3_pssm = pssm_B, pssm_weight = pssm_weight)

  df_paired <- pair_alpha_beta_multi(
    file.path(step_dir, "TRAV_TRAJ_cdr3.csv"),
    file.path(step_dir, "TRBV_TRBJ_cdr3.csv"),
    file.path(step_dir, "chainA_B_random_pair.csv")
  )
  prepare_model_csv(df_paired, peptide, mhc, file.path(step_dir, "model.csv"))

  file.path(step_dir, "model.csv")
}


# =============================================================================
# Main pipeline
# =============================================================================

run_motif_builder <- function(peptide, mhc,
                               label            = paste0(mhc, "_", peptide),
                               base_output_dir  = BASE_OUTPUT_DIR,
                               input_dir        = INPUT_DIR,
                               cdr3_baseline,
                               known_binders_file,
                               validation_file,
                               n_iterations     = N_ITERATIONS,
                               n_pairs          = N_PAIRS,
                               n_cdr3_multi     = N_CDR3_MULTI,
                               init_perc_rank   = INIT_PERC_RANK,
                               decay_factor     = DECAY_FACTOR,
                               pssm_weight      = PSSM_WEIGHT,
                               min_tcrs_pssm    = MIN_TCRS_PSSM,
                               vj_baseline_prior,
                               vj_prior_strength          = VJ_PRIOR_STRENGTH,
                               score_col                  = SCORE_COL,
                               plot_motif                 = FALSE,
                               plot_final_motif           = TRUE,
                               validate_each_step         = FALSE,
                               len_dist_conditioned_on_vj = LEN_DIST_COND_VJ) {

  step_counts <- list()

  # ------------------------------------------------------------------
  # Threshold schedule: INIT_PERC_RANK × DECAY_FACTOR^step
  # thresholds[1] = step 0, thresholds[2] = step 1, ...
  # ------------------------------------------------------------------
  thresholds <- init_perc_rank * decay_factor^(0:n_iterations)
  message(sprintf("[%s] perc_rank schedule: %s",
                  label,
                  paste(sprintf("%.3f", thresholds), collapse = " -> ")))

  # ------------------------------------------------------------------
  # Build (or reload) TEMPO predictor once
  # ------------------------------------------------------------------
  predictor_rds <- file.path(BASE_OUTPUT_DIR, label, paste0("predictor_", label, ".rds"))
  predictor_rds <- build_tempo_predictor(known_binders_file, predictor_rds)

  # ------------------------------------------------------------------
  # Step 0: flat random TCR generation — no priors
  # ------------------------------------------------------------------
  message(sprintf("\n[%s] Step 0: random TCR generation (flat baseline)", label))

  step0_dir <- file.path(base_output_dir, label, "step0")
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

  # Score step 0
  message(sprintf("[%s] Step 0: scoring (threshold = %.3f)", label, thresholds[1]))
  scored0 <- score_step(
    predictor_rds   = predictor_rds,
    model_csv       = file.path(step0_dir, "model.csv"),
    threshold       = thresholds[1],
    step            = 0,
    label           = label,
    base_output_dir = base_output_dir,
    score_col       = score_col,
    plot_motif      = plot_motif
  )
  auc01_step0 <- if (validate_each_step) {
    validate_motif(
      motif_file      = file.path(step0_dir, "model.csv"),
      validation_file = validation_file,
      output_dir      = file.path(step0_dir, "validation"),
      peptide         = peptide, mhc = mhc, plot = FALSE
    )$auc01
  } else NA_real_

  step_counts[[1]] <- data.frame(step    = 0L,
                                  n_total = scored0$n_total,
                                  n_top   = scored0$n_top,
                                  auc01   = auc01_step0)
  top_tcrs          <- scored0$top_tcrs
  all_tcrs          <- scored0$all_tcrs
  all_tcrs_combined <- scored0$all_tcrs
  current_model_csv <- file.path(step0_dir, "model.csv")
  last_step         <- 0L

  # ------------------------------------------------------------------
  # Steps 1..N_ITERATIONS: enrich → generate → score
  # ------------------------------------------------------------------
  for (iter in seq_len(n_iterations)) {

    if (nrow(top_tcrs) < 10) {
      warning(sprintf("[%s] Too few top binders after step %d — stopping early.",
                      label, iter - 1L))
      break
    }

    message(sprintf("[%s] Step %d / %d: enrichment + generation", label, iter, n_iterations))

    new_model_csv <- enrich_generate_step(
      top_tcrs                   = top_tcrs,
      all_tcrs                   = all_tcrs,
      step                       = iter,
      label                      = label,
      peptide                    = peptide, mhc = mhc,
      base_output_dir            = base_output_dir,
      cdr3_baseline              = cdr3_baseline,
      n_pairs                    = n_pairs,
      n_cdr3_multi               = n_cdr3_multi,
      pssm_weight                = pssm_weight,
      min_tcrs_pssm              = min_tcrs_pssm,
      vj_baseline_prior          = vj_baseline_prior,
      vj_prior_strength          = vj_prior_strength,
      len_dist_conditioned_on_vj = len_dist_conditioned_on_vj
    )
    if (is.null(new_model_csv)) break

    message(sprintf("[%s] Step %d / %d: scoring (threshold = %.3f)",
                    label, iter, n_iterations, thresholds[iter + 1]))

    scored_i <- score_step(
      predictor_rds   = predictor_rds,
      model_csv       = new_model_csv,
      threshold       = thresholds[iter + 1],
      step            = iter,
      label           = label,
      base_output_dir = base_output_dir,
      score_col       = score_col,
      plot_motif      = plot_motif
    )
    auc01_iter <- if (validate_each_step) {
      validate_motif(
        motif_file      = new_model_csv,
        validation_file = validation_file,
        output_dir      = file.path(base_output_dir, label,
                                     sprintf("step%d", iter), "validation"),
        peptide         = peptide, mhc = mhc, plot = FALSE
      )$auc01
    } else NA_real_

    step_counts[[iter + 1L]] <- data.frame(step    = iter,
                                            n_total = scored_i$n_total,
                                            n_top   = scored_i$n_top,
                                            auc01   = auc01_iter)
    top_tcrs          <- scored_i$top_tcrs
    all_tcrs          <- scored_i$all_tcrs
    all_tcrs_combined <- rbind(all_tcrs_combined, scored_i$all_tcrs)
    current_model_csv <- new_model_csv
    last_step         <- iter
  }

  # ------------------------------------------------------------------
  # Final motif plots (optimizer mode): top binders from last iteration
  # plotted against (1) all TCRs from that iteration and (2) standard
  # MixTCRviz baseline.
  # ------------------------------------------------------------------
  if (plot_final_motif && nrow(top_tcrs) >= 10) {
    final_step_dir <- file.path(base_output_dir, label, sprintf("step%d", last_step))
    plot_iter_motif(top_tcrs, label, last_step,
                    output_dir = file.path(final_step_dir, "motif_vs_iter"),
                    all_tcrs   = all_tcrs_combined)
    plot_iter_motif(top_tcrs, label, last_step,
                    output_dir = file.path(final_step_dir, "motif_vs_baseline"))
  }

  # ------------------------------------------------------------------
  # Final validation
  # ------------------------------------------------------------------
  message(sprintf("[%s] Final validation (AUC)", label))

  auc_result <- validate_motif(
    motif_file      = current_model_csv,
    validation_file = validation_file,
    output_dir      = file.path(base_output_dir, label, "validation"),
    peptide         = peptide, mhc = mhc,
    plot            = plot_motif
  )

  list(
    final_auc    = auc_result$auc,
    final_auc01  = auc_result$auc01,
    final_roc    = auc_result$roc,
    final_pred   = auc_result$pred,
    step_counts  = do.call(rbind, step_counts)
  )
}


### ---- Run ------------------------------------------------------------------
# Guard: only execute when run directly, not when sourced by the optimizer.
# Capture then immediately reset the flag so it never carries over between runs.

baseline      <- MixTCRviz::baseline_HomoSapiens
cdr3_baseline <- baseline$countCDR3.VJL
vj_baseline_prior <- list(A = extract_vj_baseline_prior("A"), B = extract_vj_baseline_prior("B"))

.run_pipeline         <- !exists(".sourced_by_optimizer") || !.sourced_by_optimizer
.sourced_by_optimizer <- FALSE   # reset so future direct runs always work

if (.run_pipeline) {

  auc_summary <- list()

  for (epitope in epitopes) {
    mhc     <- sub("_.*", "", epitope)
    peptide <- sub("^[^_]*_", "", epitope)
    message(sprintf("\n========== %s ==========", epitope))

    res <- run_motif_builder(
      peptide            = peptide,
      mhc                = mhc,
      cdr3_baseline      = cdr3_baseline,
      vj_baseline_prior  = vj_baseline_prior,
      known_binders_file = file.path(BASE_OUTPUT_DIR, epitope, paste0(epitope, ".csv")),
      validation_file    = file.path(BASE_OUTPUT_DIR, epitope, "validation.csv")
    )

    auc_summary[[epitope]] <- list(auc = res$final_auc, auc01 = res$final_auc01)
  }

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
