library(TEMPOtrain)
library(MixTCRviz)
library(dplyr)
library(tidyr)
library(pROC)

# =============================================================================
# ESM Motif Builder
#
# Structure-based iterative TCR motif building pipeline.
# Uses ESMFold2 ipTM scoring with dummy partner chains.
#
# Each generated batch is written twice:
#   1. model_<chain>.csv       — V/J gene names + CDR3 + metadata (as before)
#   2. model_<chain>_seqs.csv  — full amino-acid sequences for ESMFold's fold.py
#                                (ID, TCRA, TCRB, MHC, B2M, PEPTIDE); a local
#                                port of build_input.py (see Module 0).
#
# Pipeline structure:
#   Step 0:    Generate random TCRs (flat baseline) paired with dummy chain
#              → model_alpha.csv / model_alpha_seqs.csv  (real α + dummy β)
#              → model_beta.csv  / model_beta_seqs.csv   (dummy α + real β)
#              → Run ESMFold (fold.py) on cluster; place results as:
#                  step0/output_alpha.csv
#                  step0/output_beta.csv
#   Steps 1–N: Load ESM scores → filter (iptm_pair_mean >= threshold)
#              → enrich V/J + CDR3 length + PSSM per chain
#              → generate new batches → run ESMFold on cluster
#   "final":   Load last ESM scores → filter → train TEMPO on top binders
#              (dummy chain columns cleared) → AUC / AUC0.1 validation
#
# Usage: set STEP before sourcing / running this script.
#   STEP = 0              → generate step-0 TCR batches
#   STEP = 1 .. N_STEPS   → enrich from step (STEP-1) ESM scores; generate new batch
#   STEP = c(1, 2, 3)     → run multiple enrichment steps in sequence
#   STEP = "final"        → TEMPO validation using step N_STEPS ESM scores
# =============================================================================

### ---- Configuration --------------------------------------------------------

STEP            <- 0 #c(1, 2, 3)        # 0, 1, ..., N_STEPS, or "final"
N_STEPS         <- 5         # total number of enrichment steps after step 0

INPUT_DIR       <- "/Users/roessner/Documents/PostDoc/Data/MixTCRviz/data_raw/CDR123/HomoSapiens"
BASE_OUTPUT_DIR <- "TCR_motif_atlas" 
SCORE_COL       <- "iptm_pair_mean"   # column in ESMFold output.txt; higher = better

# Threshold schedule: one value per step 1..N_STEPS (TCRs with score >= threshold pass)
#ESM_THRESHOLDS  <- c(0.5, 0.6, 0.7) #c(0.5, 0.6, 0.7, 0.7)
ESM_THRESHOLDS  <- c(0.5, 0.6, 0.7)

N_PAIRS          <- 400    # V/J pairs sampled per chain from top-binder distribution
N_CDR3_MULTI     <- 3      # CDR3 sequences sampled per V/J pair (enrichment steps)
MUT_WEIGHT       <- 0.1    # IC-adjusted random-mutation rate: 0 = off; >0 mixes a uniform
                           # component into each CDR3 position, strongest at low-IC (variable)
                           # sites, ~0 at conserved anchors (affinity maturation). The PSSM
                           # blend is fixed to full weight at the junction (former PSSM_WEIGHT=1).
MIN_TCRS_PSSM    <- 30     # min top binders required to build a PSSM
VJ_PRIOR_STRENGTH  <- 30   # alpha: repertoire-prior pseudocount weight for V/J shrinkage
                           # (in evidence-count units; higher = more shrinkage toward baseline)
LEN_PRIOR_STRENGTH <- 20   # beta: same idea, for the CDR3-length enrichment shrinkage

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
# Most common V/J pair + polyGly CDR3, used as a neutral scaffold so ESMFold
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
# Module 0 — Full-sequence (fold.py) input construction
#
# Local port of build_input.py. For every generated batch we additionally write
# a CSV that ESMFold's fold.py consumes directly:
#     ID, TCRA, TCRB, MHC, B2M, PEPTIDE      (full amino-acid sequences)
# TCR chains are rebuilt as  V(before CDR3) + CDR3 + J(after CDR3)  using the
# CDR123 gene reference tables; MHC/B2M sequences come from the alleles table.
#
# Dummy V genes ("*_dummy") are handled as on the cluster
# (Script_build_struct_inputs_mod.R / get_TRseq): the base gene is looked up and
# its CDR1/CDR2 loops are masked to same-length glycine runs, so the dummy chain
# keeps a neutral framework while its (polyGly) CDR3 is still inserted.
# =============================================================================

MIXTCRVIZ_DATA <- "/Users/roessner/Documents/PostDoc/Data/MixTCRviz/data_raw/CDR123"
ALLELES_TAB    <- "/Users/roessner/Documents/PostDoc/Data/preMSA/Alleles_tab.txt"

# Allele name -> AA sequence (MHC + B2M lookups).
load_alleles_seq <- function(path = ALLELES_TAB) {
  tab <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                    quote = "", comment.char = "")
  setNames(tab$seq, tab$Allele)
}

# Load TRAV/TRAJ/TRBV/TRBJ reference tables for one species. Gaps ("g"/"-") are
# stripped from the sequence columns so CDR1/CDR2/CDR3 match within 'full'.
load_gene_tables <- function(species, data_dir = MIXTCRVIZ_DATA) {
  strip_gaps <- function(x) gsub("-", "", gsub("g", "", x))
  tabs <- list()
  for (chain in c("A", "B")) for (gt in c("V", "J")) {
    key <- paste0("TR", chain, gt)
    df  <- read.csv(file.path(data_dir, species, paste0(key, ".csv")),
                    check.names = FALSE, stringsAsFactors = FALSE)
    names(df)[1] <- "name"
    for (col in intersect(c("full", "CDR1", "CDR2", "CDR3"), names(df)))
      df[[col]] <- strip_gaps(df[[col]])
    tabs[[key]] <- df
  }
  tabs
}

.re_escape   <- function(x) gsub("([][{}()*+?.\\^$|])", "\\\\\\1", x)
# 0-based string positions (to mirror Python's rfind/find), -1 if not found.
.last_fixed  <- function(hay, needle) {
  m <- gregexpr(needle, hay, fixed = TRUE)[[1]]; if (m[1] == -1L) -1L else max(m) - 1L
}
.first_fixed <- function(hay, needle) {
  m <- regexpr(needle, hay, fixed = TRUE);       if (m == -1L)    -1L else as.integer(m) - 1L
}
.lookup      <- function(tab, key) {
  if (!is.null(key) && !is.na(key) && key %in% names(tab)) unname(tab[[key]]) else NA_character_
}

# Part of a V/J gene sequence outside the CDR3 (V: before, J: after). Mirrors
# build_input.py::get_seq_outside_cdr3, plus dummy CDR1/CDR2 glycine masking.
seq_outside_cdr3 <- function(gene_type, gene_name, cdr3, gene_table) {
  is_dummy <- grepl("_dummy$", gene_name)
  idx      <- match(sub("_dummy$", "", gene_name), gene_table$name)
  if (is.na(idx)) return(NA_character_)
  full_seq  <- gene_table$full[idx]
  gene_cdr3 <- gene_table$CDR3[idx]

  if (is_dummy && gene_type == "V") {
    for (loop_col in c("CDR1", "CDR2")) {
      loop <- gene_table[[loop_col]][idx]
      if (!is.na(loop) && nzchar(loop))
        full_seq <- sub(loop, strrep("G", nchar(loop)), full_seq, fixed = TRUE)
    }
  }

  if (gene_type == "V") {
    result <- sub(paste0(.re_escape(gene_cdr3), "$"), "", full_seq)
    if (identical(result, full_seq)) {
      for (n in 5:1) {
        pos <- .last_fixed(full_seq, substr(cdr3, 1, n))
        if (pos >= 0 && nchar(full_seq) - pos <= 15) { result <- substr(full_seq, 1, pos); break }
      }
    }
  } else {
    result <- sub(paste0("^", .re_escape(gene_cdr3)), "", full_seq)
    if (identical(result, full_seq)) {
      L <- nchar(cdr3)
      for (n in 5:1) {
        pos <- .first_fixed(full_seq, substr(cdr3, L - n + 1, L))
        if (pos >= 0 && pos + n <= 15) { result <- substr(full_seq, pos + n + 1, nchar(full_seq)); break }
      }
    }
  }
  if (identical(result, full_seq)) NA_character_ else result
}

# Full TCR chain sequence: V_before_CDR3 + CDR3 + J_after_CDR3.
build_tcr_chain <- function(row, chain, gene_tables) {
  cdr3 <- row[[paste0("cdr3_TR", chain)]]
  if (is.null(cdr3) || is.na(cdr3)) return(NA_character_)
  v_seq <- seq_outside_cdr3("V", row[[paste0("TR", chain, "V")]], cdr3, gene_tables[[paste0("TR", chain, "V")]])
  j_seq <- seq_outside_cdr3("J", row[[paste0("TR", chain, "J")]], cdr3, gene_tables[[paste0("TR", chain, "J")]])
  if (is.na(v_seq) || is.na(j_seq)) return(NA_character_)
  paste0(v_seq, cdr3, j_seq)
}

.normalize_mhc <- function(name) ifelse(grepl("^[ABC][0-9]", name), paste0("HLA_", name), name)
.infer_species <- function(mhc)  if (grepl("^H2", mhc)) "MusMusculus" else "HomoSapiens"
.infer_b2m     <- function(sp)   if (sp == "MusMusculus") "Mouse_B2M" else "B2M"

# Write the fold.py input CSV (ID, TCRA, TCRB, MHC, B2M, PEPTIDE) from a model
# table (the same data.frame written to model_<chain>.csv).
write_fold_input_csv <- function(model_df, output_file,
                                 alleles = ALLELES_SEQ, gene_tables_by_species = GENE_TABLES) {
  rows <- vector("list", nrow(model_df)); n_bad <- 0L
  for (i in seq_len(nrow(model_df))) {
    row    <- as.list(model_df[i, ])
    tcr_id <- if (!is.null(row$id)) row$id else sprintf("tcr%04d", i)
    mhc    <- .normalize_mhc(if (!is.null(row$MHC_allele_a)) row$MHC_allele_a else row$MHC)
    sp     <- if (!is.null(row$species)) row$species else .infer_species(mhc)
    gt     <- gene_tables_by_species[[sp]]
    if (is.null(gt)) { n_bad <- n_bad + 1L; next }

    tcra     <- build_tcr_chain(row, "A", gt)
    tcrb     <- build_tcr_chain(row, "B", gt)
    mhc_seq  <- .lookup(alleles, mhc)
    b2m_name <- if (!is.null(row$MHC_allele_b)) row$MHC_allele_b else .infer_b2m(sp)
    b2m_seq  <- if (tolower(b2m_name) == "none") "" else .lookup(alleles, b2m_name)
    peptide  <- if (!is.null(row$peptide)) row$peptide else row$epitope

    if (is.na(tcra) || is.na(tcrb) || is.na(mhc_seq) || is.na(b2m_seq) ||
        is.null(peptide) || is.na(peptide) || !grepl("^[A-Z]+$", paste0(tcra, tcrb))) {
      n_bad <- n_bad + 1L; next
    }
    rows[[i]] <- data.frame(ID = tcr_id, TCRA = tcra, TCRB = tcrb,
                            MHC = mhc_seq, B2M = b2m_seq, PEPTIDE = peptide,
                            stringsAsFactors = FALSE)
  }
  out_df <- do.call(rbind, rows)
  if (is.null(out_df))
    out_df <- data.frame(ID = character(), TCRA = character(), TCRB = character(),
                         MHC = character(), B2M = character(), PEPTIDE = character())
  write.csv(out_df, output_file, row.names = FALSE)
  message(sprintf("    Wrote %d fold-input complexes to %s (%d skipped)",
                  nrow(out_df), basename(output_file), n_bad))
  invisible(out_df)
}


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
  write_fold_input_csv(df_A, sub("\\.csv$", "_seqs.csv", output_file))
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
  write_fold_input_csv(df_B, sub("\\.csv$", "_seqs.csv", output_file))
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
                                    n = 5, len_dist = NULL,
                                    cdr3_pssm = NULL, mut_weight = 0) {
  key <- paste0(v_seg, "_", j_seg)
  if (!key %in% names(cdr3_baseline[[chain]])) return(character(0))

  vj_counts         <- cdr3_baseline[[chain]][[key]]
  lengths_available <- names(vj_counts)

  # Sample n lengths WITH replacement from the marginal enriched length
  # distribution (excess-over-background shrunk to the baseline length prior),
  # so a strongly enriched length yields multiple CDR3s. Weights are enrichment,
  # not baseline abundance, so rare-enriched lengths are favoured.
  active_len_dist <- len_dist
  if (!is.null(active_len_dist)) {
    valid_lens <- intersect(names(active_len_dist), lengths_available)
    if (length(valid_lens) > 0) {
      probs         <- unlist(active_len_dist[valid_lens])
      probs         <- probs / sum(probs)
      selected_lens <- sample(valid_lens, size = n, replace = TRUE, prob = probs)
    } else {
      selected_lens <- sample(lengths_available, size = n, replace = TRUE)
    }
  } else {
    selected_lens <- sample(lengths_available, size = n, replace = TRUE)
  }

  seqs <- vapply(selected_lens, function(len_name) {
    count_mat <- as.matrix(vj_counts[[len_name]])
    prob_mat  <- apply(count_mat, 2, function(col) {
      if (sum(col) == 0) rep(1 / nrow(count_mat), nrow(count_mat))
      else               col / sum(col)
    })

    # Per-position baseline information content, used by BOTH the PSSM blend and
    # the IC-adjusted mutation. Scale-free (computed from the normalized baseline
    # column): IC_baseline(pos) = log2(K) - H(pos), K = alphabet size.
    max_ic <- log2(nrow(prob_mat))                     # log2(K); K = 20 for AAs
    ic_pos <- apply(prob_mat, 2, function(col) {
      p <- col[col > 0]
      max_ic + sum(p * log2(p))                        # = max_ic - H(col)
    })

    # (3) PSSM blend at FULL weight (former PSSM_WEIGHT = 1, hardcoded): at low-IC
    # junction positions the PSSM replaces the baseline; at high-IC anchors the
    # weight -> 0, so the conserved C...F framework stays baseline-governed.
    #   w(pos) = 1 - IC_baseline(pos) / log2(K)
    if (!is.null(cdr3_pssm) && !is.null(cdr3_pssm[[len_name]])) {
      pssm      <- cdr3_pssm[[len_name]]
      aa_common <- intersect(rownames(prob_mat), rownames(pssm))
      if (length(aa_common) >= 10 && ncol(pssm) == ncol(prob_mat)) {
        w_pos <- 1 - ic_pos / max_ic                   # per-position PSSM weight
        w_pos <- pmax(0, pmin(1, w_pos))               # numerical safety
        for (pos in seq_len(ncol(prob_mat))) {
          prob_mat[aa_common, pos] <- (1 - w_pos[pos]) * prob_mat[aa_common, pos] +
                                            w_pos[pos]  * pssm[aa_common, pos]
        }
        prob_mat <- apply(prob_mat, 2, function(col) col / sum(col))
      }
    }

    # (4) IC-adjusted random mutation: mix a uniform component into each position
    # so residues ABSENT from the V/J baseline can be proposed. IC-modulated like
    # the PSSM blend (strong at the variable junction, ~0 at anchors), gated on
    # the *baseline* IC so a peaked PSSM from few top binders can't suppress it.
    # Selection (ESMFold) decides which mutations persist. mut_weight = 0 = off.
    if (mut_weight > 0) {
      unif  <- rep(1 / nrow(prob_mat), nrow(prob_mat))  # uniform over AA alphabet
      m_pos <- mut_weight * (1 - ic_pos / max_ic)       # per-position mutation rate
      m_pos <- pmax(0, pmin(mut_weight, m_pos))
      for (pos in seq_len(ncol(prob_mat))) {
        prob_mat[, pos] <- (1 - m_pos[pos]) * prob_mat[, pos] + m_pos[pos] * unif
      }
      prob_mat <- apply(prob_mat, 2, function(col) col / sum(col))
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
                                     n = 5, len_dist = NULL,
                                     cdr3_pssm = NULL, mut_weight = 0) {
  df           <- read.csv(pair_file)
  chain_letter <- sub("^TR", "", chain)
  v_col        <- paste0("TR", chain_letter, "V")
  j_col        <- paste0("TR", chain_letter, "J")

  results <- mapply(draw_random_cdr3_multi,
                    chain = chain, v_seg = df[[v_col]], j_seg = df[[j_col]],
                    MoreArgs = list(cdr3_baseline = cdr3_baseline,
                                   n           = n,
                                   len_dist    = len_dist,
                                   cdr3_pssm   = cdr3_pssm,
                                   mut_weight  = mut_weight),
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
# Module 4 — ESM score loading and filtering
# =============================================================================

load_esm_scores <- function(scores_file, model_file, score_col = SCORE_COL) {
  scores <- read.csv(scores_file, stringsAsFactors = FALSE, check.names = FALSE)
  model  <- read.csv(model_file,  stringsAsFactors = FALSE, check.names = FALSE)

  # The ESMFold output (output_<chain>.csv) keys on 'ID'; the model table
  # (model_<chain>.csv) keys on 'id'. Normalise both to lowercase 'id'.
  fix_id <- function(df) {
    hit <- which(tolower(names(df)) == "id")
    if (length(hit) >= 1) names(df)[hit[1]] <- "id"
    df
  }
  scores <- fix_id(scores)
  model  <- fix_id(model)

  if (!score_col %in% colnames(scores))
    stop("Column '", score_col, "' not found in ESM scores file.\n",
         "Available: ", paste(colnames(scores), collapse = ", "))
  if (!"id" %in% colnames(scores) || !"id" %in% colnames(model))
    stop("Both the ESM scores file and model table must have an 'id'/'ID' column.")

  merge(model, scores[, c("id", score_col)], by = "id")
}


# =============================================================================
# Module 5 — Enrichment helpers (per-chain)
# =============================================================================

# Functional (gene_type == "F") gene names per gene, from the same reference
# tables step 0 uses. Non-functional genes (pseudogenes / ORF / orphons) are
# present in the baseline but must never enter the enrichment steps.
load_functional_genes <- function(input_dir) {
  genes <- c("TRAV", "TRAJ", "TRBV", "TRBJ")
  setNames(lapply(genes, function(g) {
    df <- read.csv(file.path(input_dir, paste0(g, ".csv")))
    df[df$gene_type == "F", 1]
  }), genes)
}

# Marginal repertoire priors P_baseline(V) and P_baseline(J) per chain, from the
# joint countVJ matrix (row/col sums), restricted to functional genes.
extract_vj_marginal_prior <- function(chain_letter, functional) {
  m  <- baseline_HomoSapiens$countVJ[[paste0("TR", chain_letter)]]
  pv <- rowSums(m); pj <- colSums(m)
  fv <- functional[[paste0("TR", chain_letter, "V")]]
  fj <- functional[[paste0("TR", chain_letter, "J")]]
  pv <- pv[names(pv) %in% fv & pv > 0]
  pj <- pj[names(pj) %in% fj & pj > 0]
  list(V = pv / sum(pv), J = pj / sum(pj))
}

# Marginal baseline CDR3-length prior P_baseline(L) per chain, pooled over all
# V/J pairs in countCDR3.VJL. Keys are "L_<n>" (matching lengths_available).
extract_len_baseline_prior <- function(chain_letter) {
  vj_list <- baseline_HomoSapiens$countCDR3.VJL[[paste0("TR", chain_letter)]]
  len_counts <- list()
  for (key in names(vj_list)) for (l in names(vj_list[[key]])) {
    cnt <- sum(vj_list[[key]][[l]][, 1])   # n CDR3s of length l for this V/J
    len_counts[[l]] <- (if (is.null(len_counts[[l]])) 0 else len_counts[[l]]) + cnt
  }
  v <- unlist(len_counts); v <- v[v > 0]
  as.list(v / sum(v))
}


# Factorized, evidence-weighted, baseline-shrunk V/J distribution.
# The joint (V,J) enrichment is empirically ~fully explained by the product of
# the marginal V and J effects (interaction ≈ 1), so V and J are modelled
# INDEPENDENTLY. This credits enrichment to the gene that earned it and stops an
# enriched V from dragging an unspecific "passenger" J along (and vice versa).
# Per gene:
#   excess(g)    = max(0, n_top(g) - n_all(g) * p_global)   [excess over background]
#   posterior(g) ∝ alpha * P_baseline(g) + excess(g)        [shrink to marginal prior]
# Restricted to functional genes (marginal_prior names). Returns list(V, J) or NULL.
extract_vj_marginal_posterior <- function(top_tcrs, scored, chain_letter, marginal_prior, alpha) {
  v_col <- paste0("TR", chain_letter, "V")
  j_col <- paste0("TR", chain_letter, "J")
  if (!all(c(v_col, j_col) %in% colnames(top_tcrs))) return(NULL)
  if (nrow(scored) == 0) return(NULL)
  p_global <- nrow(top_tcrs) / nrow(scored)

  post <- function(gene_col, prior) {
    n_top <- table(top_tcrs[[gene_col]][!is.na(top_tcrs[[gene_col]])])
    n_all <- table(scored[[gene_col]][!is.na(scored[[gene_col]])])
    genes <- names(prior)
    excess <- vapply(genes, function(g) {
      nt <- if (g %in% names(n_top)) n_top[[g]] else 0
      na <- if (g %in% names(n_all)) n_all[[g]] else 0
      max(0, nt - na * p_global)
    }, numeric(1))
    posterior <- alpha * unlist(prior) + excess
    if (sum(posterior) == 0) return(NULL)
    posterior / sum(posterior)
  }

  pv <- post(v_col, marginal_prior$V)
  pj <- post(j_col, marginal_prior$J)
  if (is.null(pv) || is.null(pj)) return(NULL)
  list(V = pv, J = pj)
}


# Sample n_pairs (V, J) from the FACTORIZED distribution:
#   P_sample(v, j) ∝ posterior_V(v) * posterior_J(j) * 1[(v,j) has a CDR3 baseline]
# V and J drawn by their own marginal enrichment, combined only into pairs the
# reference repertoire actually observed (generable and, since posteriors are
# functional-only, functional).
sample_vj_pairs <- function(vj_dist, chain_letter, n_pairs, output_dir, cdr3_baseline) {
  v_col <- paste0("TR", chain_letter, "V")
  j_col <- paste0("TR", chain_letter, "J")
  if (is.null(vj_dist))
    stop(sprintf("No V/J distribution available for chain TR%s.", chain_letter))

  keys    <- names(cdr3_baseline[[paste0("TR", chain_letter)]])
  split_k <- strsplit(keys, "_")
  vg      <- vapply(split_k, `[`, character(1), 1)
  jg      <- vapply(split_k, `[`, character(1), 2)
  ok      <- vg %in% names(vj_dist$V) & jg %in% names(vj_dist$J)
  if (!any(ok))
    stop(sprintf("No valid functional V/J pairs with a CDR3 baseline for chain TR%s.", chain_letter))
  vg <- vg[ok]; jg <- jg[ok]

  w   <- as.numeric(vj_dist$V[vg]) * as.numeric(vj_dist$J[jg])   # factorized weight per pair
  idx <- sample(seq_along(w), size = n_pairs, replace = TRUE, prob = w / sum(w))

  pairs <- data.frame(vg[idx], jg[idx], stringsAsFactors = FALSE)
  colnames(pairs) <- c(v_col, j_col)
  write.csv(pairs, file.path(output_dir, paste0(v_col, "_", j_col, ".csv")),
            row.names = FALSE)
  pairs
}


# Marginal CDR3-length distribution, same excess-over-background + baseline
# shrinkage as the V/J marginals (a raw freq_top/freq_all ratio blows up at rare
# lengths and, with with-replacement length sampling, hijacks generation toward
# spurious extreme lengths):
#   posterior(L) ∝ alpha * P_baseline_len(L) + max(0, n_top(L) - n_all(L) * p_global)
# alpha = len_prior_strength. Returns a named list over "L_<n>" keys.
extract_cdr3_len_dist <- function(top_tcrs, scored, chain_letter, len_baseline_prior, alpha) {
  cdr3_col <- paste0("cdr3_TR", chain_letter)
  if (!cdr3_col %in% colnames(top_tcrs)) return(NULL)

  Lt <- nchar(top_tcrs[[cdr3_col]][!is.na(top_tcrs[[cdr3_col]])])
  La <- nchar(scored[[cdr3_col]][!is.na(scored[[cdr3_col]])])
  if (length(Lt) == 0 || length(La) == 0) return(NULL)

  n_top    <- table(paste0("L_", Lt))
  n_all    <- table(paste0("L_", La))
  p_global <- length(Lt) / length(La)

  lens   <- names(len_baseline_prior)
  excess <- vapply(lens, function(l) {
    nt <- if (l %in% names(n_top)) n_top[[l]] else 0
    na <- if (l %in% names(n_all)) n_all[[l]] else 0
    max(0, nt - na * p_global)
  }, numeric(1))
  posterior <- alpha * unlist(len_baseline_prior) + excess
  if (sum(posterior) == 0) return(NULL)
  as.list(posterior / sum(posterior))
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
  message(sprintf("[%s] → Run ESMFold on cluster, then place results as:\n  %s\n  %s",
                  label,
                  file.path(step_dir, "output_alpha.csv"),
                  file.path(step_dir, "output_beta.csv")))

  invisible(step_dir)
}


# Enrich one chain from previous ESM scores and generate a new model CSV.
# Returns top_tcrs (filtered) or NULL if too few binders.
enrich_one_chain <- function(chain_letter, step, label, peptide, mhc_allele, species,
                              base_output_dir, cdr3_baseline,
                              n_pairs, n_cdr3_multi, mut_weight, min_tcrs_pssm,
                              vj_baseline_prior, vj_prior_strength,
                              len_baseline_prior, len_prior_strength,
                              threshold) {

  chain_name    <- if (chain_letter == "A") "alpha" else "beta"
  prev_step_dir <- file.path(base_output_dir, label, sprintf("step%d", step - 1))
  step_dir      <- file.path(base_output_dir, label, sprintf("step%d", step))
  dir.create(step_dir, showWarnings = FALSE, recursive = TRUE)

  scores_file <- file.path(prev_step_dir, sprintf("output_%s.csv", chain_name))
  model_file  <- file.path(prev_step_dir, sprintf("model_%s.csv",  chain_name))

  if (!file.exists(scores_file))
    stop(sprintf("[%s] ESM scores file missing: %s\nRun ESMFold (fold.py) on the cluster first.", label, scores_file))
  if (!file.exists(model_file))
    stop(sprintf("[%s] Model table missing: %s\nThis is the model_%s.csv batch generated for this step.", label, model_file, chain_name))

  scored   <- load_esm_scores(scores_file, model_file)
  top_tcrs <- scored[scored[[SCORE_COL]] >= threshold, ]
  message(sprintf("  [chain %s] %d / %d TCRs with %s >= %.2f",
                  chain_letter, nrow(top_tcrs), nrow(scored), SCORE_COL, threshold))

  write.csv(top_tcrs, file.path(prev_step_dir, sprintf("top_tcrs_%s.csv", chain_name)),
            row.names = FALSE)

  if (nrow(top_tcrs) < 10) {
    warning(sprintf("[%s] Too few top binders for chain %s (n=%d) — skipping enrichment.",
                    label, chain_letter, nrow(top_tcrs)))
    return(NULL)
  }

  # (1) V/J distribution: factorized marginal enrichment (V and J independent),
  # each shrunk toward the functional-only marginal prior; pairs formed only from
  # combinations that have a CDR3 baseline.
  vj_dist <- extract_vj_marginal_posterior(top_tcrs, scored, chain_letter,
                                           vj_baseline_prior[[chain_letter]], vj_prior_strength)
  sample_vj_pairs(vj_dist, chain_letter, n_pairs, step_dir, cdr3_baseline)

  # (2) Marginal CDR3 length distribution (excess-over-background shrunk to the
  # baseline length prior, strength = len_prior_strength)
  len_dist <- extract_cdr3_len_dist(top_tcrs, scored, chain_letter,
                                    len_baseline_prior[[chain_letter]], len_prior_strength)

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
                           n = n_cdr3_multi, len_dist = len_dist,
                           cdr3_pssm = cdr3_pssm, mut_weight = mut_weight)

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
                             n_pairs, n_cdr3_multi, mut_weight, min_tcrs_pssm,
                             vj_baseline_prior, vj_prior_strength,
                             len_baseline_prior, len_prior_strength,
                             esm_thresholds, species = "HomoSapiens",
                             plot_each_step     = PLOT_EACH_STEP,
                             validate_each_step = VALIDATE_EACH_STEP,
                             validation_file    = NULL, mhc = NULL) {

  threshold     <- esm_thresholds[step]
  prev_step_dir <- file.path(base_output_dir, label, sprintf("step%d", step - 1))
  step_dir      <- file.path(base_output_dir, label, sprintf("step%d", step))

  message(sprintf("[%s] Step %d: enrichment + generation (threshold = %.2f)", label, step, threshold))

  top_alpha <- enrich_one_chain("A", step, label, peptide, mhc_allele, species,
                                base_output_dir, cdr3_baseline,
                                n_pairs, n_cdr3_multi, mut_weight, min_tcrs_pssm,
                                vj_baseline_prior, vj_prior_strength,
                                len_baseline_prior, len_prior_strength,
                                threshold)
  top_beta  <- enrich_one_chain("B", step, label, peptide, mhc_allele, species,
                                base_output_dir, cdr3_baseline,
                                n_pairs, n_cdr3_multi, mut_weight, min_tcrs_pssm,
                                vj_baseline_prior, vj_prior_strength,
                                len_baseline_prior, len_prior_strength,
                                threshold)

  # Optional: MixTCRviz motif plots — saved into prev_step_dir (scores source)
  if (plot_each_step) {
    if (!is.null(top_alpha) && nrow(top_alpha) >= 10)
      plot_step_motif(top_alpha, "A", label, step - 1, prev_step_dir)
    if (!is.null(top_beta) && nrow(top_beta) >= 10)
      plot_step_motif(top_beta, "B", label, step - 1, prev_step_dir)
  }

  # Always write the combined top-binder motif file into prev_step_dir
  model_label <- paste0(mhc, "_", peptide)
  motif_df    <- build_motif_df(top_alpha, top_beta, model_label)
  motif_file  <- NULL
  if (!is.null(motif_df)) {
    motif_file <- file.path(prev_step_dir, "motif_top_binders.csv")
    write.csv(motif_df, motif_file, row.names = FALSE)
  }

  # Optional: TEMPO validation — saved into prev_step_dir (scores source)
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
  message(sprintf("[%s] → Run ESMFold on cluster, then place results as:\n  %s\n  %s",
                  label,
                  file.path(step_dir, "output_alpha.csv"),
                  file.path(step_dir, "output_beta.csv")))

  invisible(list(top_alpha = top_alpha, top_beta = top_beta))
}


run_final_validation <- function(label, peptide, mhc, mhc_allele,
                                  base_output_dir, n_steps, esm_thresholds,
                                  validation_file) {
  last_step_dir <- file.path(base_output_dir, label, sprintf("step%d", n_steps))
  threshold     <- esm_thresholds[n_steps]

  message(sprintf("[%s] Final: loading step-%d ESM scores (threshold = %.2f)",
                  label, n_steps, threshold))

  collect_top <- function(chain_letter) {
    chain_name  <- if (chain_letter == "A") "alpha" else "beta"
    scores_file <- file.path(last_step_dir, sprintf("output_%s.csv", chain_name))
    model_file  <- file.path(last_step_dir, sprintf("model_%s.csv",  chain_name))
    if (!file.exists(scores_file)) {
      warning(sprintf("[%s] ESM scores file missing: %s — skipping chain %s.",
                      label, scores_file, chain_letter))
      return(NULL)
    }
    scored   <- load_esm_scores(scores_file, model_file)
    top_tcrs <- scored[scored[[SCORE_COL]] >= threshold, ]
    message(sprintf("  [chain %s] %d / %d top binders", chain_letter, nrow(top_tcrs), nrow(scored)))
    write.csv(top_tcrs, file.path(last_step_dir, sprintf("top_tcrs_%s.csv", chain_name)),
              row.names = FALSE)
    top_tcrs
  }

  top_alpha <- collect_top("A")
  top_beta  <- collect_top("B")

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

baseline         <- MixTCRviz::baseline_HomoSapiens
cdr3_baseline    <- baseline$countCDR3.VJL
functional_genes <- load_functional_genes(INPUT_DIR)   # gene_type == "F" only
vj_baseline_prior  <- list(A = extract_vj_marginal_prior("A", functional_genes),
                           B = extract_vj_marginal_prior("B", functional_genes))
len_baseline_prior <- list(A = extract_len_baseline_prior("A"),
                           B = extract_len_baseline_prior("B"))

# Reference data for the full-sequence (fold.py) input CSVs (Module 0).
ALLELES_SEQ <- load_alleles_seq()
GENE_TABLES <- list()
for (sp in c("HomoSapiens", "MusMusculus")) {
  if (dir.exists(file.path(MIXTCRVIZ_DATA, sp))) GENE_TABLES[[sp]] <- load_gene_tables(sp)
}

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
        mut_weight                 = MUT_WEIGHT,
        min_tcrs_pssm              = MIN_TCRS_PSSM,
        vj_baseline_prior          = vj_baseline_prior,
        vj_prior_strength          = VJ_PRIOR_STRENGTH,
        len_baseline_prior         = len_baseline_prior,
        len_prior_strength         = LEN_PRIOR_STRENGTH,
        esm_thresholds             = ESM_THRESHOLDS,
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
      esm_thresholds  = ESM_THRESHOLDS,
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
