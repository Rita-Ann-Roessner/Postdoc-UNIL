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

STEP            <- 1        # 0, 1, ..., N_STEPS, or "final"
N_STEPS         <- 2         # total number of enrichment steps after step 0

INPUT_DIR       <- "/Users/roessner/Documents/PostDoc/Data/MixTCRviz/data_raw/CDR123/HomoSapiens"
BASE_OUTPUT_DIR <- "TCR_motif_atlas_cdr3_templated_dana_orig_VJ" 
SCORE_COL       <- "iptm_pair_mean"   # column in ESMFold output.txt; higher = better

# Threshold schedule: one value per step 1..N_STEPS (TCRs with score >= threshold pass)
ESM_THRESHOLDS  <- c(0.5, 0.6)

N_PAIRS          <- 400    # V/J pairs sampled per chain from top-binder distribution
N_CDR3_MULTI     <- 3      # CDR3 sequences sampled per V/J pair (enrichment steps)
MIN_TCRS_PSSM    <- 30     # min top binders required to build a PSSM
VJ_PRIOR_STRENGTH <- 60    # α: repertoire-prior pseudocount weight for the V/J posterior,
                           # in units of TCR counts -- roughly how much evidence a gene needs
                           # before it moves off the repertoire baseline. The value used for
                           # TCR_motif_atlas_cdr3_templated.
                           # A multiplicative rate-ratio form was tried and removed: at the gene
                           # level it was too noisy, because with ~47 TRAV / 61 TRAJ / 66 TRBV
                           # genes and a few hundred passers many genes sit at low n_all, where a
                           # ratio is unstable. Additive keeps rare genes pinned near the prior,
                           # since their absolute excess is necessarily small. LENGTH still uses
                           # the ratio form (LEN_WEIGHT_K) -- there the clipping was the binding
                           # problem, not the noise.
LEN_WEIGHT_K     <- 5      # pseudocount in the per-length enrichment weight
                           #   w(L) = (n_top(L) + k) / (n_all(L)*p_global + k)
                           # Pulls lengths with few observations toward neutral (w=1).
                           # Replaces the former LEN_PRIOR_STRENGTH / LEN_ENRICH_MAX_FOLD:
                           # the ratio is intrinsically bounded, so no cap is needed.

# --- Germline-templated CDR3 generation (Module 3b) --------------------------
# TRUE  = build each CDR3 as  V-germline block + junction + J-germline block,
#         with the trimming extents sampled from the V/J/L baseline PWM itself.
#         Germline blocks are emitted ATOMICALLY, so e.g. TRAJ30's "DDKIIF"
#         can never be mangled into "DGKIIF" by independent per-position draws.
# FALSE = legacy behaviour (every position drawn independently).
CDR3_TEMPLATED  <- TRUE
JUNCTION_PSSM   <- FALSE   # TRUE = let the top-binder PSSM shape the junction residues.
                           # FALSE = junction drawn from the V/J baseline only, so no
                           # top-binder information enters the CDR3 sequence at all
                           # (enrichment then acts solely via V/J and length choice).
TEMPLATE_MAX_TRIES <- 20   # redraws allowed when v_ext + j_ext would exceed L

# Where the trimming extents come from (Module 3a):
#   "survival" = extent distribution read off the V/J/L baseline PWM alone.
#   "figure"   = same distribution, but SHIFTED so its mode lands on the modal
#                post-trimming germline extent for that chain and length, as
#                measured off the iReceptor panel of the V/J-vs-CDR3-length
#                figure (FIG_EXT below). Spread and per-gene conditioning are
#                preserved; only the centre moves.
# Measured effect of "figure" over the step-0 panel:
#   alpha  V+J 80.7% -> 87.8% of the CDR3 (junction 2.68 -> 1.69 residues)
#   beta   V+J 65.9% -> 71.4%             (junction 5.08 -> 4.26 residues)
# J moves ~2.7x as much as V on both chains, which is the intended bias.
EXTENT_SOURCE <- "figure"

# --- Fixed step-0 panel ------------------------------------------------------
# Path to a one-off generated step-0 repertoire reused for EVERY epitope, so that
# tcrNNNN is literally the same TCR everywhere and step-0 scores are directly
# paired across peptides (differences are peptide effects, not sampling noise).
# Only peptide / MHC / species are re-stamped per epitope; the TCRs are identical.
# Create it once with:  Rscript generate_step0_panel.R
# NULL = draw a fresh step-0 repertoire per epitope (the original behaviour).
STEP0_PANEL <- "step0_panel"

# Decoy V/J background correction (optional). When enabled, the null background in
# the V/J excess is each gene's pass propensity against a DECOY (non-cognate)
# epitope — the peptide-INDEPENDENT ESMFold fold artifact — instead of the flat
# global pass rate. Passing above the decoy level counts as enrichment; germline
# fold artifacts (e.g. TRAV12-2) stay at baseline, while genuinely epitope-specific
# V/J (high cognate, low decoy) are still enriched. Recomputed on the fly at each
# step's own selection threshold from the decoy step-0 fold.
# Every gene present in the decoy fold uses its own rate, however few observations
# back it; only genes absent from the decoy fall back to p_global.
DECOY_CORRECTION_VJ  <- FALSE                       # correct the V/J excess null background
DECOY_DIR        <- "step0_background"             # holds <decoy_label>/step0 folds
DECOY_BY_MHC     <- c(A0201 = "A0201_ALAAAAAAV")   # MHC (allele w/o HLA_ prefix) -> decoy label

FINAL_PERCENTILE <- 95  # final (STEP="final") top-binder cutoff = this percentile of the POOLED
                        # (alpha+beta) score distribution of the LAST step — i.e. keep the top
                        # (100 - N)% best-scoring TCRs of the final step. Epitope-adaptive.
                        # NULL = use the fixed ESM_THRESHOLDS[N_STEPS] instead.

CHAIN_WEIGHT <- "linear" # how TEMPO combines the two chains when scoring validation:
                        # "equal"  = default TEMPO score (score_A + score_B, equal weight).
                        # "sqrt"   = weight each chain by sqrt(n_passers): score = sqrt(n_A)*score_A
                        #            + sqrt(n_B)*score_B (gentle; trusts the better-sampled chain).
                        # "linear" = weight by the passer count itself: n_A*score_A + n_B*score_B
                        #            (aggressive; can over-penalise a real but minority chain).

PLOT_EACH_STEP     <- TRUE   # if TRUE, plot MixTCRviz motif after each enrichment step
VALIDATE_EACH_STEP <- TRUE  # if TRUE, run TEMPO validation after each enrichment step

VALIDATION_AVAILABLE <- TRUE  # set FALSE when this epitope has NO validation data — skips ALL
                              # TEMPO validation (per-step AND final); motifs are still built and
                              # written. Manual flag on purpose (not auto file-detection), so a
                              # misplaced validation.csv can't silently pass unnoticed.

EPITOPES_FILE <- NULL #"epitopes.txt"   # one "MHC_PEPTIDE" label per line, e.g. "A0201_ELAGIGILTV"
epitopes <- c("A0201_ELAGIGILTV", "A0201_GILGFVFTL", "A0201_LLWNGPMAV")

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

  # germline-templated assembly (same path the enrichment steps use), so step 0
  # and steps 1+ produce CDR3s of the same construction. Falls through to plain
  # per-position sampling when no template is available.
  if (isTRUE(CDR3_TEMPLATED)) {
    s <- assemble_templated_cdr3(chain, v_seg, j_seg, prob_mat)
    if (!is.na(s)) return(s)
  }

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

# `weights` (one per sequence, default 1) lets the decoy correction reach the PSSM:
# a CDR3 from a gene that passes at its decoy rate contributes little or nothing.
# All-1 weights reproduce the uncorrected counts exactly.
build_cdr3_pssm <- function(cdr3_seqs, weights = NULL, pseudocount = 0.1) {
  if (is.null(weights)) weights <- rep(1, length(cdr3_seqs))
  keep      <- !is.na(cdr3_seqs) & nchar(cdr3_seqs) > 0
  cdr3_seqs <- cdr3_seqs[keep]; weights <- weights[keep]
  if (length(cdr3_seqs) == 0) return(NULL)

  amino_acids <- c("A","C","D","E","F","G","H","I","K","L",
                   "M","N","P","Q","R","S","T","V","W","Y")
  idx_by_length <- split(seq_along(cdr3_seqs), nchar(cdr3_seqs))

  pssm_list <- lapply(names(idx_by_length), function(len_str) {
    idx <- idx_by_length[[len_str]]
    L   <- as.integer(len_str)
    mat <- matrix(pseudocount, nrow = length(amino_acids), ncol = L,
                  dimnames = list(amino_acids, NULL))
    for (i in idx) {
      chars <- strsplit(cdr3_seqs[i], "")[[1]]
      for (pos in seq_along(chars)) {
        aa <- chars[pos]
        if (aa %in% amino_acids) mat[aa, pos] <- mat[aa, pos] + weights[i]
      }
    }
    apply(mat, 2, function(col) col / sum(col))
  })
  names(pssm_list) <- paste0("L_", names(idx_by_length))
  pssm_list
}


# =============================================================================
# Module 3a — Germline CDR3 templates and trimming extents
#
# A rearranged CDR3 is  V-germline(trimmed) + N-additions(+D) + J-germline(trimmed).
# The germline halves are deterministic given the gene; what varies is HOW MUCH of
# each survives exonuclease trimming. That trimming profile is already recoverable
# from the baseline PWM: at a position still covered by the V template, the PWM
# shows the germline residue with high probability, and that probability decays as
# you walk inward. So
#       P(V extent >= i)  ~=  PWM[ v_germline[i], i ]
# is a survival curve, and differencing it gives the extent distribution. Same from
# the other end for J. No external recombination model (OLGA/IGoR) needed.
# =============================================================================

# gene -> its germline CDR3 contribution ("CAVN" for TRAV12-2, "NRDDKIIF" for TRAJ30)
load_germline_cdr3 <- function(input_dir = INPUT_DIR) {
  strip_gaps <- function(x) gsub("-", "", gsub("g", "", x))
  setNames(lapply(c("TRAV", "TRAJ", "TRBV", "TRBJ"), function(seg) {
    df <- read.csv(file.path(input_dir, paste0(seg, ".csv")), stringsAsFactors = FALSE)
    names(df)[1] <- "gene"
    v <- strip_gaps(df$CDR3)
    setNames(v, df$gene)[nzchar(v) & !is.na(v)]
  }), c("TRAV", "TRAJ", "TRBV", "TRBJ"))
}

# Template lookup with composite-name fallback: "TRBV6-2/6-3" is a merge of two
# genes that are identical over CDR1/2/3, and only the composite name is absent
# from the reference table, so resolving to the part before "/" loses nothing.
germline_cdr3_part <- function(seg, gene, germline = GERMLINE_CDR3) {
  tab <- germline[[seg]]
  if (is.null(tab) || is.null(gene) || is.na(gene)) return(NA_character_)
  if (gene %in% names(tab)) return(unname(tab[[gene]]))
  if (grepl("/", gene, fixed = TRUE)) {
    base <- sub("/.*$", "", gene)
    if (base %in% names(tab)) return(unname(tab[[base]]))
  }
  NA_character_
}

# Modal post-trimming germline extents by chain and CDR3 length, read off the
# iReceptor panel (A) of the V/J-contribution-vs-CDR3-length figure by measuring
# each bubble and taking the row-wise argmax. These are POST-trimming counts, so
# they are directly comparable to v_ext / j_ext.
#
# Confidence is uneven. L_Ja is the only panel with real signal (R = 0.446) and
# iReceptor and VDJdb agree on it residue-for-residue across all nine lengths;
# L_Va agrees between datasets to within one residue. The two beta columns rest
# on iReceptor alone -- VDJdb's beta panels are saturated (all three columns
# measure 16-17 px) and cannot be read. Beta's J values also sit above what the
# validation positives show, possibly an annotation-convention difference in how
# ambiguous junction residues are split between D and J; revisit if beta motifs
# degrade. Outside L = 10-18 the bubbles are too small to measure, so L is
# clamped rather than extrapolated.
FIG_EXT <- list(
  A = list(V = c(`10`=3,`11`=3,`12`=3,`13`=3,`14`=4,`15`=5,`16`=5,`17`=5,`18`=5),
           J = c(`10`=7,`11`=7,`12`=8,`13`=9,`14`=9,`15`=10,`16`=11,`17`=12,`18`=12)),
  B = list(V = c(`10`=4,`11`=4,`12`=4,`13`=4,`14`=4,`15`=5,`16`=6,`17`=6,`18`=6),
           J = c(`10`=6,`11`=6,`12`=6,`13`=6,`14`=6,`15`=7,`16`=7,`17`=8,`18`=8))
)

# Target (v_ext, j_ext) for one chain letter and CDR3 length; NULL if unknown.
fig_extent <- function(chain_letter, L) {
  tb <- FIG_EXT[[chain_letter]]
  if (is.null(tb)) return(NULL)
  Lc <- as.character(max(10L, min(18L, as.integer(L))))
  c(v = unname(tb$V[[Lc]]), j = unname(tb$J[[Lc]]))
}

# P(extent = i) for i = 1..length(surv), from a survival curve P(extent >= i).
# cummin enforces the nesting constraint (if position 3 is templated, 1 and 2 are
# too) and incidentally denoises the PWM estimate. Extent 0 is excluded on
# purpose: the conserved C (V-encoded) and F/W (J-encoded) are the definitional
# boundaries of a CDR3 and are never trimmed away, so every CDR3 keeps at least
# one germline residue at each end.
.extent_probs <- function(surv) {
  surv <- cummin(pmax(0, pmin(1, surv)))
  p <- surv - c(surv[-1], 0)
  if (sum(p) <= 0) return(NULL)
  p / sum(p)
}

# Sample how far the V and J germline templates survive, given the baseline PWM
# for this V/J pair at this length. Returns NULL if no valid (non-overlapping)
# pair of extents can be drawn.
sample_extents <- function(prob_mat, v_tmpl, j_tmpl, L,
                           max_tries = TEMPLATE_MAX_TRIES,
                           chain_letter = NULL) {
  aas <- rownames(prob_mat)
  surv_at <- function(res, pos) if (res %in% aas && pos >= 1 && pos <= L) prob_mat[res, pos] else 0

  # Cap each extent at L-1 so neither template can consume the whole sequence:
  # guarantees (v_ext = 1, j_ext = 1) is always available and the retry loop below
  # terminates. L >= 2 is required for both anchors to exist at all.
  if (L < 2) return(NULL)
  vc <- strsplit(v_tmpl, "")[[1]]; jc <- strsplit(j_tmpl, "")[[1]]
  nv <- min(length(vc), L - 1);    nj <- min(length(jc), L - 1)
  if (nv < 1 || nj < 1) return(NULL)

  surv_v <- vapply(seq_len(nv), function(i) surv_at(vc[i], i), numeric(1))
  # J is aligned to the C-terminus: its last residue sits at position L
  surv_j <- vapply(seq_len(nj), function(k)
                     surv_at(jc[length(jc) - k + 1], L - k + 1), numeric(1))

  pv <- .extent_probs(surv_v); pj <- .extent_probs(surv_j)
  if (is.null(pv) || is.null(pj)) return(NULL)

  # EXTENT_SOURCE == "figure": offset each draw so the distribution's mode lands
  # on the figure's modal extent for this chain/length. The PWM still supplies
  # the shape (and with it the per-gene trimming spread); only the centre moves.
  # A target the gene cannot reach -- Jalpha = 12 against a 9-residue TRAJ12
  # template -- is absorbed by the clamp to nv / nj.
  off_v <- 0L; off_j <- 0L
  if (identical(EXTENT_SOURCE, "figure") && !is.null(chain_letter)) {
    tgt <- fig_extent(chain_letter, L)
    if (!is.null(tgt)) {
      off_v <- as.integer(tgt[["v"]]) - which.max(pv)
      off_j <- as.integer(tgt[["j"]]) - which.max(pj)
    }
  }
  draw_v <- function() max(1L, min(nv, sample(seq_len(nv), 1, prob = pv) + off_v))
  draw_j <- function() max(1L, min(nj, sample(seq_len(nj), 1, prob = pj) + off_j))

  for (i in seq_len(max_tries)) {
    v_ext <- draw_v()
    j_ext <- draw_j()
    if (v_ext + j_ext <= L) return(list(v_ext = v_ext, j_ext = j_ext))
  }
  # persistent overlap (short L, long templates): keep V, clamp J to what fits,
  # never below 1 so the C-terminal F/W anchor stays germline
  v_ext <- min(draw_v(), L - 1)
  list(v_ext = v_ext, j_ext = max(1, min(nj, L - v_ext)))
}

# Assemble ONE CDR3 as  V-germline block + junction + J-germline block, given the
# baseline PWM for this V/J pair at this length. The germline blocks are copied
# verbatim (so they are emitted atomically); only the junction is sampled per
# position. Used by BOTH step 0 (draw_random_cdr3) and the enrichment steps
# (draw_random_cdr3_multi) so the two stay consistent.
#   base_prob  baseline PWM -> trimming survival curves
#   junc_src   PWM the junction residues are drawn from (= base_prob unless the
#              caller wants a PSSM-blended junction)
# Returns NA_character_ when no template exists or the extents are degenerate;
# the caller then falls back to plain per-position sampling.
assemble_templated_cdr3 <- function(chain, v_seg, j_seg, base_prob, junc_src = base_prob) {
  cl     <- sub("^TR", "", chain)
  v_tmpl <- germline_cdr3_part(paste0("TR", cl, "V"), v_seg)
  j_tmpl <- germline_cdr3_part(paste0("TR", cl, "J"), j_seg)
  if (is.na(v_tmpl) || is.na(j_tmpl)) return(NA_character_)

  L   <- ncol(base_prob)
  ext <- sample_extents(base_prob, v_tmpl, j_tmpl, L, chain_letter = cl)
  if (is.null(ext)) return(NA_character_)
  n_junc <- L - ext$v_ext - ext$j_ext
  if (n_junc < 0) return(NA_character_)

  junc <- character(0)
  if (n_junc > 0) {
    aas  <- rownames(junc_src)
    junc <- vapply(seq_len(n_junc), function(k) {
      p <- junc_src[, ext$v_ext + k]
      if (sum(p) == 0) NA_character_ else sample(aas, 1, prob = p)
    }, character(1))
    if (any(is.na(junc))) return(NA_character_)
  }
  paste0(substr(v_tmpl, 1, ext$v_ext),
         paste0(junc, collapse = ""),
         substring(j_tmpl, nchar(j_tmpl) - ext$j_ext + 1))
}


# =============================================================================
# Module 3 — Multi-CDR3 sampling with V/J + length + PSSM enrichments
# =============================================================================

draw_random_cdr3_multi <- function(chain, v_seg, j_seg, cdr3_baseline,
                                    n = 5, len_weight = NULL,
                                    cdr3_pssm = NULL, mut_weight = 0) {
  key <- paste0(v_seg, "_", j_seg)
  if (!key %in% names(cdr3_baseline[[chain]])) return(character(0))

  vj_counts         <- cdr3_baseline[[chain]][[key]]
  lengths_available <- names(vj_counts)

  # ---- length: THIS pair's own baseline, tilted by the per-length enrichment ---
  # CDR3 length is largely determined by the J gene (per-J mean CDR3a length spans
  # 10.1-14.0 while within-J sd is 0.4-1.6), so a length distribution pooled over
  # all V/J pairs is badly mis-conditioned for any single pair. Use the pair's own
  # baseline (exactly what step 0 does) and tilt it by how much better that length
  # performed than the pool average (extract_cdr3_len_weight):
  #     P(L | V,J)  ∝  count_pair(L)  *  w(L)
  # w = 1 is neutral, so with no enrichment this reduces to step 0's behaviour.
  # Lengths absent from the weight table get 1 rather than being dropped.
  pair_cnt <- vapply(lengths_available,
                     function(l) sum(vj_counts[[l]][, 1], na.rm = TRUE), numeric(1))
  names(pair_cnt) <- lengths_available
  valid_lens <- lengths_available[pair_cnt > 0]

  if (length(valid_lens) > 0) {
    w <- rep(1, length(valid_lens)); names(w) <- valid_lens
    if (!is.null(len_weight)) {
      for (l in valid_lens)
        if (l %in% names(len_weight)) w[[l]] <- as.numeric(len_weight[[l]])
    }
    probs <- pair_cnt[valid_lens] * w
  } else {
    probs <- numeric(0)
  }

  if (length(valid_lens) > 0 && sum(probs) > 0) {
    probs         <- probs / sum(probs)
    # n DISTINCT lengths (no replacement). This caps any single length at 1/n of
    # the batch, so a strongly concentrated preference cannot be reproduced at
    # generation time -- but it also damps the step-to-step feedback loop.
    #
    # Sampling WITH replacement was tried (TCR_motif_atlas_cdr3_templated_mod_length)
    # and reverted: the realised distribution then follows `probs` faithfully, but
    # `probs` is built from the PREVIOUS step's pass rates, so generation chases a
    # one-step-stale, noisy ranking. For LLW alpha that meant 264 -> 45 TCRs at
    # L=10 (the reference mode, and the best-performing length at 8.3%) and
    # 189 -> 529 at the worst lengths, cutting the alpha pass rate 4.6% -> 2.2%.
    # Keeping replacement off trades peak concentration for stability.
    selected_lens <- sample(valid_lens, size = min(n, sum(probs > 0)),
                            replace = FALSE, prob = probs)
  } else {
    # pair has no usable length counts at all -> uniform over whatever it supports
    selected_lens <- sample(lengths_available, size = min(n, length(lengths_available)))
  }

  seqs <- vapply(selected_lens, function(len_name) {
    count_mat <- as.matrix(vj_counts[[len_name]])
    prob_mat  <- apply(count_mat, 2, function(col) {
      if (sum(col) == 0) rep(1 / nrow(count_mat), nrow(count_mat))
      else               col / sum(col)
    })
    # untouched baseline copy: the trimming survival curves must be read off the
    # BASELINE, not off a PSSM-blended matrix (and it is also what the junction is
    # drawn from when JUNCTION_PSSM = FALSE).
    base_prob <- prob_mat

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

    # (4) IC-adjusted Dirichlet perturbation: randomly re-emphasise amino acids
    # ALREADY present in the blend by drawing the per-position distribution from
    # Dirichlet(alpha * P_blend). Support is preserved (an AA with p = 0 gets
    # shape 0 -> stays 0), so this diversifies WITHIN the observed alphabet only —
    # no novel residues, no uniform noise. IC-gated on the *baseline* IC like the
    # PSSM blend (strong at the variable junction, ~0 at conserved anchors) so a
    # peaked PSSM from few top binders can't over-concentrate the junction.
    # Selection (ESMFold) decides which variants persist. mut_weight = 0 = off.
    #   s_pos = mut_weight * (1 - IC/log2K)     per-position perturbation strength
    #   alpha = (1 - s_pos) / s_pos             Dirichlet concentration:
    #             s->0  => alpha->Inf => distribution unchanged (anchors)
    #             s=0.1 => alpha=9    => mild reshuffle around the blend
    #             s->1  => alpha->0   => single blend-weighted residue (max diversify)
    if (mut_weight > 0) {
      s_pos <- mut_weight * (1 - ic_pos / max_ic)
      s_pos <- pmax(0, pmin(1, s_pos))
      for (pos in seq_len(ncol(prob_mat))) {
        s <- s_pos[pos]
        if (s <= 0) next                                    # anchor: no perturbation
        p     <- prob_mat[, pos]
        alpha <- (1 - s) / s                                # total Dirichlet concentration
        g     <- rgamma(length(p), shape = alpha * p, rate = 1)  # p_i = 0 -> 0 (support kept)
        if (sum(g) > 0) {
          prob_mat[, pos] <- g / sum(g)
        } else {
          # alpha -> 0 numerical underflow: realise the Dirichlet limit explicitly
          # (all mass on one residue, drawn proportional to the blend).
          pick <- sample.int(length(p), 1, prob = p)
          v <- numeric(length(p)); v[pick] <- 1
          prob_mat[, pos] <- v
        }
      }
      prob_mat <- apply(prob_mat, 2, function(col) col / sum(col))
    }

    amino_acids <- rownames(prob_mat)
    L <- ncol(prob_mat)

    # ---- (5) germline-templated assembly ------------------------------------
    # V block + junction + J block, via the same helper step 0 uses. Extents are
    # always read off the untouched baseline; only the junction may be drawn from
    # the PSSM-blended matrix, and only when JUNCTION_PSSM = TRUE.
    if (isTRUE(CDR3_TEMPLATED)) {
      s <- assemble_templated_cdr3(chain, v_seg, j_seg, base_prob,
                                   junc_src = if (isTRUE(JUNCTION_PSSM)) prob_mat else base_prob)
      if (!is.na(s)) return(s)
      # no template / degenerate extents -> fall through to per-position sampling
    }

    seq_vec <- sapply(seq_len(L), function(pos) {
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
                                     n = 5, len_weight = NULL,
                                     cdr3_pssm = NULL, mut_weight = 0) {
  df           <- read.csv(pair_file)
  chain_letter <- sub("^TR", "", chain)
  v_col        <- paste0("TR", chain_letter, "V")
  j_col        <- paste0("TR", chain_letter, "J")

  results <- mapply(draw_random_cdr3_multi,
                    chain = chain, v_seg = df[[v_col]], j_seg = df[[j_col]],
                    MoreArgs = list(cdr3_baseline = cdr3_baseline,
                                   n          = n,
                                   len_weight = len_weight,
                                   cdr3_pssm  = cdr3_pssm,
                                   mut_weight = mut_weight),
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


# Metadata table that matches a folded output_<chain>.csv. Prefer the IMMUTABLE
# input_<chain>.csv (the exact batch that was sent to ESMFold) when present, so
# scores are never joined against a model_<chain>.csv that a later run may have
# regenerated (positional tcrNNNN ids -> silent score/TCR mismatch). Falls back
# to model_<chain>.csv when no frozen input exists.
meta_csv <- function(step_dir, chain_name) {
  inp <- file.path(step_dir, sprintf("input_%s.csv", chain_name))
  if (file.exists(inp)) inp else file.path(step_dir, sprintf("model_%s.csv", chain_name))
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


# Per-V and per-J DECOY pass propensity for one chain, pooled over the decoy step-0
# replicates (DECOY_DIR/rep*/<decoy_label>/step0), computed at `threshold`. The
# decoy is a non-cognate peptide, so this is the peptide-INDEPENDENT ESMFold fold
# artifact per gene. `mhc` selects the decoy via DECOY_BY_MHC (HLA_ prefix
# stripped). Returns list(V, J = pass-rate vectors; n_V, n_J = sample counts), or
# NULL if no decoy is configured/folded for this MHC.
# Step-0 fold directories holding one background label. Supports both layouts:
#   <dir>/<label>/step0          single fold of the fixed step-0 panel (current)
#   <dir>/rep<N>/<label>/step0   independent CDR3 replicates (legacy)
# Returns every directory that exists, so a mixed tree still works -- though old
# replicate folds predate germline templating and should not be pooled with new ones.
background_step0_dirs <- function(dir, label) {
  flat <- file.path(dir, label, "step0")
  reps <- list.dirs(dir, recursive = FALSE)
  reps <- reps[grepl("^rep[0-9]+$", basename(reps))]
  cand <- c(flat, file.path(reps, label, "step0"))
  cand[dir.exists(cand)]
}

load_decoy_background <- function(chain_letter, mhc, threshold,
                                  decoy_dir = DECOY_DIR, decoy_by_mhc = DECOY_BY_MHC,
                                  score_col = SCORE_COL) {
  mhc_key <- sub("^HLA_", "", mhc)
  if (!mhc_key %in% names(decoy_by_mhc)) return(NULL)
  label      <- decoy_by_mhc[[mhc_key]]
  chain_name <- if (chain_letter == "A") "alpha" else "beta"
  v_col      <- paste0("TR", chain_letter, "V")
  j_col      <- paste0("TR", chain_letter, "J")

  frames <- list()
  for (sd in background_step0_dirs(decoy_dir, label)) {
    sf <- file.path(sd, sprintf("output_%s.csv", chain_name))
    mf <- file.path(sd, sprintf("model_%s.csv",  chain_name))
    if (file.exists(sf) && file.exists(mf))
      frames[[length(frames) + 1]] <- load_esm_scores(sf, mf, score_col)[, c(v_col, j_col, score_col)]
  }
  if (length(frames) == 0) return(NULL)

  d    <- do.call(rbind, frames)                # pool all replicate decoy TCRs
  pass <- d[[score_col]] >= threshold
  rate_n <- function(key) {
    keep <- !is.na(key)
    a <- table(key[keep]); p <- table(key[keep & pass]); g <- names(a)
    list(rate = setNames(as.numeric(ifelse(g %in% names(p), p[g], 0)) / as.numeric(a[g]), g),
         n    = setNames(as.integer(a[g]), g))
  }
  V <- rate_n(d[[v_col]]); J <- rate_n(d[[j_col]])
  list(V = V$rate, n_V = V$n, J = J$rate, n_J = J$n)
}


# Factorized, evidence-weighted, baseline-shrunk V/J distribution.
# The joint (V,J) enrichment is empirically ~fully explained by the product of
# the marginal V and J effects (interaction ≈ 1), so V and J are modelled
# INDEPENDENTLY. This credits enrichment to the gene that earned it and stops an
# enriched V from dragging an unspecific "passenger" J along (and vice versa).
# Per gene the null background is p_global (pooled pass rate) by default, OR -- when
# a decoy background is supplied (decoy_bg from load_decoy_background) -- that gene's
# own decoy pass propensity, so a peptide-independent fold artifact isn't credited
# as enrichment:
#   bg(g)        = decoy_bg(g)  if supplied, else p_global
#   excess(g)    = max(0, n_top(g) - n_all(g)*bg(g))     [passers above the null]
#   posterior(g) ∝ α*P_baseline(g) + excess(g)           [shrink toward the prior]
#
# Known costs of this form, measured at GIL step 1 -- kept deliberately, but do not
# re-derive them from scratch:
#  - BELOW background max(0,.) erases performance entirely (39 of 45 alpha genes),
#    leaving those genes ranked purely by naive repertoire frequency: TRAV13-1
#    passed 0/18 and got 2.8% of the posterior while TRAV27 -- GIL's correct gene --
#    passed 1/9 (at background) and got 1.2%.
#  - ABOVE background `excess` is a COUNT, so it scales with the slots a gene already
#    holds: posterior ∝ slots * rate_advantage, i.e. geometric growth. TRAV12-2 ran
#    51 -> 261 -> 525 slots while TRAV27 collapsed to 9.
#  - α is a total pseudocount spread over the prior, so α/n_genes gives TRBJ
#    (13 genes) ~3.6x more shrinkage per gene than TRAV (45 genes).
# A rate-ratio form fixed all three but was too noisy at the gene level; see
# VJ_PRIOR_STRENGTH. Length still uses the ratio (extract_cdr3_len_weight).
# decoy_bg = NULL falls back to p_global for every gene.
# Restricted to functional genes (marginal_prior names). Returns list(V, J) or NULL.
extract_vj_marginal_posterior <- function(top_tcrs, scored, chain_letter, marginal_prior,
                                           alpha = VJ_PRIOR_STRENGTH, decoy_bg = NULL) {
  v_col <- paste0("TR", chain_letter, "V")
  j_col <- paste0("TR", chain_letter, "J")
  if (!all(c(v_col, j_col) %in% colnames(top_tcrs))) return(NULL)
  if (nrow(scored) == 0) return(NULL)
  p_global <- nrow(top_tcrs) / nrow(scored)

  post <- function(gene_col, prior, bg_p) {
    n_top <- table(top_tcrs[[gene_col]][!is.na(top_tcrs[[gene_col]])])
    n_all <- table(scored[[gene_col]][!is.na(scored[[gene_col]])])
    genes <- names(prior)
    # counts and the per-gene null, shared by both modes
    nt_v <- vapply(genes, function(g) if (g %in% names(n_top)) as.numeric(n_top[[g]]) else 0, numeric(1))
    na_v <- vapply(genes, function(g) if (g %in% names(n_all)) as.numeric(n_all[[g]]) else 0, numeric(1))
    # decoy pass propensity as the null; fall back to the pooled global rate
    # only when no decoy background is supplied for this gene.
    bg_v <- vapply(genes, function(g)
      if (!is.null(bg_p) && g %in% names(bg_p)) as.numeric(bg_p[[g]]) else p_global, numeric(1))

    excess    <- pmax(0, nt_v - na_v * bg_v)
    posterior <- alpha * unlist(prior) + excess
    if (sum(posterior) == 0) return(NULL)
    posterior / sum(posterior)
  }

  pv <- post(v_col, marginal_prior$V, decoy_bg$V)
  pj <- post(j_col, marginal_prior$J, decoy_bg$J)
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


# Per-length enrichment WEIGHT: how much better or worse CDR3s of length L passed
# than the pool average. Multiplies a V/J pair's own baseline length distribution
# in draw_random_cdr3_multi, so w = 1 means "this length performed exactly at the
# pool average -> leave the pair's prior alone".
#
#   w(L) = (n_top(L) + k) / (n_all(L) * p_global + k)
#
# This replaces the former excess-over-background posterior
#   max(0, n_top - n_all*p_global), normalised against P_baseline_len,
# whose max(0, .) clipped every below-average length to exactly 0. Those lengths
# then all collapsed to one identical weight (0.35 for alpha, 0.31 for beta),
# discarding the ordering across ~half the data -- and always the LONG half, since
# the pass rate falls monotonically with length. That asymmetry is what drove the
# ~1.5-residue-per-step drift toward short CDR3s.
#
# The ratio keeps the ordering on both sides of the average, is intrinsically
# bounded (no cap needed), and needs no repertoire baseline: the pair's own
# counts already supply that in draw_random_cdr3_multi.
# k is a pseudocount so a length with few observations is pulled toward neutral
# (alpha L=9 with 2/5 would otherwise read 3.5x; with k=5 it reads 1.26x).
# `weights` (one per top-binder row, default 1) carries the decoy correction in:
# passers from decoy-prone genes count for less. p_global is weighted the same
# way so the comparison stays on a consistent scale. All-1 weights = uncorrected.
# Returns a named list over "L_<n>" keys.
extract_cdr3_len_weight <- function(top_tcrs, scored, chain_letter, k = LEN_WEIGHT_K,
                                    weights = NULL) {
  cdr3_col <- paste0("cdr3_TR", chain_letter)
  if (!cdr3_col %in% colnames(top_tcrs)) return(NULL)
  if (is.null(weights)) weights <- rep(1, nrow(top_tcrs))

  keep_t <- !is.na(top_tcrs[[cdr3_col]])
  Lt <- nchar(top_tcrs[[cdr3_col]][keep_t]); wt <- weights[keep_t]
  La <- nchar(scored[[cdr3_col]][!is.na(scored[[cdr3_col]])])
  if (length(Lt) == 0 || length(La) == 0) return(NULL)

  n_top    <- tapply(wt, paste0("L_", Lt), sum)
  n_all    <- table(paste0("L_", La))
  p_global <- sum(wt) / length(La)

  lens <- union(names(n_all), names(n_top))
  w <- vapply(lens, function(l) {
    nt <- if (l %in% names(n_top)) as.numeric(n_top[[l]]) else 0
    na <- if (l %in% names(n_all)) as.numeric(n_all[[l]]) else 0
    (nt + k) / (na * p_global + k)
  }, numeric(1))
  if (!length(w)) return(NULL)
  as.list(w)
}


# =============================================================================
# Module 6 — TEMPO final validation
# Trains TEMPO directly from the top-binder motif file (no percentile-rank
# predictor needed); uses raw score for ROC.
# =============================================================================

validate_motif <- function(motif_file, validation_file, output_dir,
                            peptide, mhc, score_col = "score", plot = TRUE,
                            chain_weight = CHAIN_WEIGHT) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Score the validation set with a given chain model; use TEMPO's returned data
  # frame directly (res$score with Label + score) — no disk glob, so there is no
  # stale-file / locale-sort ambiguity.
  score_chain <- function(ch) {
    res <- TEMPOtrain(
      input.train = motif_file, input.pred = validation_file, output.path = output_dir,
      chain = ch, build.prank = FALSE, compute.prank = FALSE,
      write.data.pred = FALSE, write.data.train = FALSE, write.predictor = FALSE,
      return.object = TRUE)
    res$score
  }

  if (chain_weight %in% c("sqrt", "linear")) {
    # Weight chains by their motif passer counts: score = wA*score_A + wB*score_B,
    # with wA/wB = n (linear) or sqrt(n) (gentler).
    motif <- read.csv(motif_file)
    nA <- sum(!is.na(motif$cdr3_TRA)); nB <- sum(!is.na(motif$cdr3_TRB))
    wfun <- if (chain_weight == "sqrt") sqrt else identity
    wA <- wfun(nA); wB <- wfun(nB)
    dA <- score_chain("A"); dB <- score_chain("B")
    if (nrow(dA) != nrow(dB))
      stop("Per-chain prediction row counts differ (A=", nrow(dA), ", B=", nrow(dB), ").")
    pred_df <- dA
    pred_df[[score_col]] <- wA * dA[[score_col]] + wB * dB[[score_col]]
    message(sprintf("  chain weights (%s): A=%.2f (n=%d)  B=%.2f (n=%d)", chain_weight, wA, nA, wB, nB))
  } else {
    # equal weight = default combined (AB) TEMPO score
    pred_df <- score_chain("AB")
  }

  # Persist the final scored validation set (with the combined/weighted score in
  # score_col) as a single file — written from memory, so no glob ambiguity.
  write.csv(pred_df, file.path(output_dir, "validation_pred.csv"), row.names = FALSE)

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
                       base_output_dir, input_dir, species = "HomoSapiens",
                       panel_dir = STEP0_PANEL) {
  step_dir <- file.path(base_output_dir, label, "step0")
  dir.create(step_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- fixed panel: reuse one repertoire for every epitope -------------------
  # The TCRs (and their tcrNNNN ids) are identical across epitopes; only the
  # peptide / MHC / species columns are re-stamped. Nothing is re-drawn.
  if (!is.null(panel_dir) && nzchar(panel_dir)) {
    message(sprintf("[%s] Step 0: using fixed panel '%s'", label, panel_dir))
    for (cn in c("alpha", "beta")) {
      src <- file.path(panel_dir, sprintf("model_%s.csv", cn))
      if (!file.exists(src))
        stop(sprintf("STEP0_PANEL is set but '%s' is missing.\nCreate the panel once with: Rscript generate_step0_panel.R", src))
      df <- read.csv(src, stringsAsFactors = FALSE)
      df$peptide <- peptide; df$MHC <- mhc_allele; df$species <- species
      out <- file.path(step_dir, sprintf("model_%s.csv", cn))
      write.csv(df, out, row.names = FALSE)
      write_fold_input_csv(df, sub("\\.csv$", "_seqs.csv", out))
    }
    n_alpha <- nrow(read.csv(file.path(step_dir, "model_alpha.csv")))
    n_beta  <- nrow(read.csv(file.path(step_dir, "model_beta.csv")))
    message(sprintf("[%s] Step 0 done: %d alpha-batch TCRs, %d beta-batch TCRs (from panel).",
                    label, n_alpha, n_beta))
    message(sprintf("[%s] → Run ESMFold on cluster, then place results as:\n  %s\n  %s",
                    label, file.path(step_dir, "output_alpha.csv"),
                    file.path(step_dir, "output_beta.csv")))
    return(invisible(step_dir))
  }

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
                              n_pairs, n_cdr3_multi, min_tcrs_pssm,
                              vj_baseline_prior,
                              threshold, mut_weight = 0) {

  chain_name    <- if (chain_letter == "A") "alpha" else "beta"
  prev_step_dir <- file.path(base_output_dir, label, sprintf("step%d", step - 1))
  step_dir      <- file.path(base_output_dir, label, sprintf("step%d", step))
  dir.create(step_dir, showWarnings = FALSE, recursive = TRUE)

  scores_file <- file.path(prev_step_dir, sprintf("output_%s.csv", chain_name))
  model_file  <- meta_csv(prev_step_dir, chain_name)   # prefer frozen input_<chain>.csv

  if (!file.exists(scores_file))
    stop(sprintf("[%s] ESM scores file missing: %s\nRun ESMFold (fold.py) on the cluster first.", label, scores_file))
  if (!file.exists(model_file))
    stop(sprintf("[%s] Metadata table missing: %s\nExpected the frozen input_%s.csv (or model_%s.csv) matching the folded output.", label, model_file, chain_name, chain_name))

  scored   <- load_esm_scores(scores_file, model_file)
  top_tcrs <- scored[scored[[SCORE_COL]] >= threshold, ]
  message(sprintf("  [chain %s] %d / %d TCRs with %s >= %.2f",
                  chain_letter, nrow(top_tcrs), nrow(scored), SCORE_COL, threshold))

  write.csv(top_tcrs, file.path(prev_step_dir, sprintf("top_tcrs_%s.csv", chain_name)),
            row.names = FALSE)
  # Also write the top binders in fold.py sequence format (ID, TCRA, TCRB, MHC,
  # B2M, PEPTIDE) so they can be re-folded directly, same as model_<chain>_seqs.csv.
  write_fold_input_csv(top_tcrs, file.path(prev_step_dir, sprintf("top_tcrs_%s_seqs.csv", chain_name)))

  if (nrow(top_tcrs) < 10) {
    warning(sprintf("[%s] Too few top binders for chain %s (n=%d) — skipping enrichment.",
                    label, chain_letter, nrow(top_tcrs)))
    return(NULL)
  }

  # Decoy V/J background at THIS step's threshold (peptide-independent fold
  # artifact), loaded when DECOY_CORRECTION_VJ is on.
  decoy_bg <- NULL
  if (isTRUE(DECOY_CORRECTION_VJ)) {
    decoy_bg <- load_decoy_background(chain_letter, mhc_allele, threshold)
    if (is.null(decoy_bg))
      warning(sprintf("[%s] DECOY_CORRECTION_VJ on but no decoy background found for MHC %s (chain %s) — using p_global.",
                      label, mhc_allele, chain_letter))
  }

  # (1) V/J distribution: factorized marginal enrichment (V and J independent),
  # each shrunk toward the functional-only marginal prior; pairs formed only from
  # combinations that have a CDR3 baseline. When DECOY_CORRECTION_VJ is on, the
  # excess null background is each gene's decoy pass propensity instead of p_global.
  vj_dist <- extract_vj_marginal_posterior(top_tcrs, scored, chain_letter,
                                           vj_baseline_prior[[chain_letter]],
                                           alpha = VJ_PRIOR_STRENGTH, decoy_bg = decoy_bg)
  sample_vj_pairs(vj_dist, chain_letter, n_pairs, step_dir, cdr3_baseline)

  # (1b) Carry the SAME correction into the length posterior and the PSSM.
  # Without this, the decoy correction only reshapes which V/J get sampled while
  # (2) and (3) are still estimated from the raw, artifact-contaminated top set.
  # Per TCR the weight is the fraction of its gene's passes that are excess over
  # that gene's decoy background -- the same quantity (1) uses, applied per row:
  #     w(g) = max(0, n_top(g) - n_all(g)*bg(g)) / n_top(g)      in [0, 1]
  # A gene passing exactly at its decoy rate contributes nothing; a gene with no
  # decoy propensity keeps full weight. V and J weights multiply, matching the
  # factorized treatment in (1). No decoy -> all weights 1 -> exact no-op.
  tcr_weights <- rep(1, nrow(top_tcrs))
  if (!is.null(decoy_bg)) {
    p_global <- nrow(top_tcrs) / nrow(scored)
    gene_w <- function(gt, bg_p) {
      col   <- paste0("TR", chain_letter, gt)
      n_all <- table(scored[[col]][!is.na(scored[[col]])])
      n_top <- table(top_tcrs[[col]][!is.na(top_tcrs[[col]])])
      vapply(top_tcrs[[col]], function(g) {
        if (is.na(g) || !(g %in% names(n_top))) return(1)
        nt <- as.numeric(n_top[[g]])
        na <- if (g %in% names(n_all)) as.numeric(n_all[[g]]) else 0
        bg <- if (!is.null(bg_p) && g %in% names(bg_p)) bg_p[[g]] else p_global
        if (nt <= 0) 0 else min(1, max(0, (nt - na * bg) / nt))
      }, numeric(1))
    }
    tcr_weights <- gene_w("V", decoy_bg$V) * gene_w("J", decoy_bg$J)
    message(sprintf("  decoy weights (chain %s): mean=%.2f  zero-weight TCRs=%d/%d",
                    chain_letter, mean(tcr_weights), sum(tcr_weights == 0), length(tcr_weights)))
  }

  # (2) Per-length enrichment weight (ratio to the pool average; see
  # extract_cdr3_len_weight). Tilts each V/J pair's own length distribution.
  len_weight <- extract_cdr3_len_weight(top_tcrs, scored, chain_letter,
                                        k = LEN_WEIGHT_K, weights = tcr_weights)

  # (3) PSSM from top binders' CDR3s
  cdr3_col <- paste0("cdr3_TR", chain_letter)
  if (nrow(top_tcrs) >= min_tcrs_pssm) {
    cdr3_pssm <- build_cdr3_pssm(top_tcrs[[cdr3_col]], weights = tcr_weights)
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
                           n = n_cdr3_multi, len_weight = len_weight,
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
                             n_pairs, n_cdr3_multi, min_tcrs_pssm,
                             vj_baseline_prior,
                             esm_thresholds, species = "HomoSapiens", mut_weight = 0,
                             plot_each_step     = PLOT_EACH_STEP,
                             validate_each_step = VALIDATE_EACH_STEP,
                             validation_file    = NULL, mhc = NULL) {

  threshold     <- esm_thresholds[step]
  prev_step_dir <- file.path(base_output_dir, label, sprintf("step%d", step - 1))
  step_dir      <- file.path(base_output_dir, label, sprintf("step%d", step))

  message(sprintf("[%s] Step %d: enrichment + generation (threshold = %.2f)", label, step, threshold))

  top_alpha <- enrich_one_chain("A", step, label, peptide, mhc_allele, species,
                                base_output_dir, cdr3_baseline,
                                n_pairs, n_cdr3_multi, min_tcrs_pssm,
                                vj_baseline_prior,
                                threshold, mut_weight)
  top_beta  <- enrich_one_chain("B", step, label, peptide, mhc_allele, species,
                                base_output_dir, cdr3_baseline,
                                n_pairs, n_cdr3_multi, min_tcrs_pssm,
                                vj_baseline_prior,
                                threshold, mut_weight)

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
  if (isTRUE(VALIDATION_AVAILABLE) && validate_each_step &&
      !is.null(motif_file) && !is.null(validation_file) && !is.null(mhc)) {
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
                                  validation_file,
                                  final_percentile = FINAL_PERCENTILE) {
  last_step_dir <- file.path(base_output_dir, label, sprintf("step%d", n_steps))

  # Final cutoff: either the fixed schedule value, or the Nth percentile of the
  # POOLED (alpha+beta) score distribution of the LAST step — keep the top
  # (100 - N)% best-scoring TCRs. Epitope-adaptive to the final pool.
  if (!is.null(final_percentile)) {
    last_scores <- unlist(lapply(c("alpha", "beta"), function(ch) {
      f <- file.path(last_step_dir, sprintf("output_%s.csv", ch))
      if (file.exists(f)) {
        d <- read.csv(f, check.names = FALSE)
        if (SCORE_COL %in% names(d)) d[[SCORE_COL]] else numeric(0)
      } else numeric(0)
    }))
    if (length(last_scores) == 0)
      stop(sprintf("[%s] Cannot compute percentile: no step-%d output scores.", label, n_steps))
    threshold <- as.numeric(quantile(last_scores, final_percentile / 100, na.rm = TRUE))
    message(sprintf("[%s] Final: threshold = %.3f (%gth pctile of pooled step-%d scores, n=%d)",
                    label, threshold, final_percentile, n_steps, length(last_scores)))
  } else {
    threshold <- esm_thresholds[n_steps]
    message(sprintf("[%s] Final: threshold = %.2f (fixed schedule)", label, threshold))
  }

  collect_top <- function(chain_letter) {
    chain_name  <- if (chain_letter == "A") "alpha" else "beta"
    scores_file <- file.path(last_step_dir, sprintf("output_%s.csv", chain_name))
    model_file  <- meta_csv(last_step_dir, chain_name)   # prefer frozen input_<chain>.csv
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
    # Also write fold.py sequence format for potential re-folding.
    write_fold_input_csv(top_tcrs, file.path(last_step_dir, sprintf("top_tcrs_%s_seqs.csv", chain_name)))
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

  if (!isTRUE(VALIDATION_AVAILABLE)) {
    message(sprintf("[%s] VALIDATION_AVAILABLE = FALSE — motif written; skipping TEMPO validation.", label))
    return(NULL)
  }

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
GERMLINE_CDR3    <- load_germline_cdr3(INPUT_DIR)      # gene -> germline CDR3 block
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

# Guard: skip the cluster-in-the-loop main run when sourced by the benchmark
# harness (which only needs the functions + setup above).
.run_main <- !exists(".sourced_for_benchmark") || !.sourced_for_benchmark
if (.run_main) {

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
        min_tcrs_pssm              = MIN_TCRS_PSSM,
        vj_baseline_prior          = vj_baseline_prior,
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

}  # end .run_main guard (skipped when sourced with .sourced_for_benchmark <- TRUE)
