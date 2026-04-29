library(MixTCRviz)

### Dummy chain definitions ###
# Most common V/J pair for TRB: TRBV20-1 + TRBJ2-7, most common CDR3 length = 14
# CDR3 format: C + 12xG + F
DUMMY_TRB <- list(
  v   = "TRBV20-1_dummy",
  j   = "TRBJ2-7",
  cdr3 = paste0("C", paste(rep("G", 12), collapse = ""), "F")
)

# Most common V/J pair for TRA: TRAV38-2/DV8 + TRAJ43, most common CDR3 length = 12
# CDR3 format: C + 10xG + F
DUMMY_TRA <- list(
  v    = "TRAV38-2/DV8_dummy",
  j    = "TRAJ43",
  cdr3 = paste0("C", paste(rep("G", 10), collapse = ""), "F")
)


### Helper Functions ###

draw_random_cdr3 <- function(chain, v_seg, j_seg, cdr3_baseline) {
  if (!chain %in% c("TRA", "TRB")) stop("chain must be 'TRA' or 'TRB'")

  key <- paste0(v_seg, "_", j_seg)
  if (!key %in% names(cdr3_baseline[[chain]])) {
    return(NA)
  }

  vj_counts <- cdr3_baseline[[chain]][[key]]
  lengths_available <- names(vj_counts)

  lengths_cnt <- sapply(lengths_available, function(len_name) {
    sum(vj_counts[[len_name]][, 1], na.rm = TRUE)
  })

  if (sum(lengths_cnt) == 0) return(NA)

  fractions <- lengths_cnt / sum(lengths_cnt)
  random_length <- sample(lengths_available, size = 1, prob = fractions)

  count_mat <- as.matrix(vj_counts[[random_length]])

  prob_mat <- apply(count_mat, 2, function(col) {
    if (sum(col) == 0) {
      rep(0, length(col))
    } else {
      col / sum(col)
    }
  })

  amino_acids <- rownames(prob_mat)
  L <- ncol(prob_mat)

  seq_vec <- sapply(1:L, function(pos) {
    probs <- prob_mat[, pos]
    if (sum(probs) == 0) return(NA)
    sample(amino_acids, 1, prob = probs)
  })

  if (any(is.na(seq_vec))) return(NA)

  paste0(seq_vec, collapse = "")
}


generate_vj_pairs <- function(chain, input_dir, output_dir) {
  genes <- c(paste0("TR", chain, "V"), paste0("TR", chain, "J"))
  lst <- vector("list", length(genes))

  for (i in seq_along(genes)) {
    gene <- genes[i]
    df_gene <- read.csv(file.path(input_dir, paste0(gene, ".csv")))

    df_gene <- df_gene[df_gene$gene_type == "F", , drop = FALSE]
    lst[[i]] <- df_gene[[1]]
  }

  pairs <- expand.grid(lst[[1]], lst[[2]])
  colnames(pairs) <- genes

  outfile <- file.path(output_dir, paste0(genes[1], "_", genes[2], ".csv"))
  write.csv(pairs, outfile, row.names = FALSE)

  pairs
}


sample_chain_cdr3 <- function(chain, pair_file, cdr3_baseline, output_file) {
  df <- read.csv(pair_file)

  if (chain == "TRA") {
    v_col <- "TRAV"
    j_col <- "TRAJ"
  } else if (chain == "TRB") {
    v_col <- "TRBV"
    j_col <- "TRBJ"
  } else {
    stop("chain must be 'TRA' or 'TRB'")
  }

  df$random_CDR3 <- mapply(
    draw_random_cdr3,
    chain = chain,
    v_seg = df[[v_col]],
    j_seg = df[[j_col]],
    MoreArgs = list(cdr3_baseline = cdr3_baseline)
  )

  df <- na.omit(df)
  write.csv(df, output_file, row.names = FALSE)

  df
}


# Pair each real alpha chain row with the fixed dummy beta chain.
pair_alpha_with_dummy_beta <- function(df_A, output_file) {
  df_A <- df_A
  colnames(df_A)[colnames(df_A) == "random_CDR3"] <- "cdr3_TRA"

  df_A$TRBV    <- DUMMY_TRB$v
  df_A$TRBJ    <- DUMMY_TRB$j
  df_A$cdr3_TRB <- DUMMY_TRB$cdr3

  write.csv(df_A, output_file, row.names = FALSE)
  df_A
}


# Pair each real beta chain row with the fixed dummy alpha chain.
pair_beta_with_dummy_alpha <- function(df_B, output_file) {
  df_B <- df_B
  colnames(df_B)[colnames(df_B) == "random_CDR3"] <- "cdr3_TRB"

  df_B$TRAV    <- DUMMY_TRA$v
  df_B$TRAJ    <- DUMMY_TRA$j
  df_B$cdr3_TRA <- DUMMY_TRA$cdr3

  write.csv(df_B, output_file, row.names = FALSE)
  df_B
}


prepare_epitope_input <- function(df_paired, peptide, mhc, species = "HomoSapiens", output_file) {
  df <- df_paired
  df$peptide  <- peptide
  df$MHC      <- mhc
  df$species  <- species

  write.csv(df, output_file, row.names = FALSE)
  df
}


run_pipeline_for_peptide <- function(peptide, mhc, base_output_dir, input_dir, cdr3_baseline) {
  topdir <- file.path(base_output_dir, paste0("step1_dummy_", peptide))
  dir.create(topdir, showWarnings = FALSE, recursive = TRUE)

  generate_vj_pairs("A", input_dir, topdir)
  generate_vj_pairs("B", input_dir, topdir)

  df_A <- sample_chain_cdr3(
    chain      = "TRA",
    pair_file  = file.path(topdir, "TRAV_TRAJ.csv"),
    cdr3_baseline = cdr3_baseline,
    output_file = file.path(topdir, "TRAV_TRAJ_cdr3.csv")
  )

  df_B <- sample_chain_cdr3(
    chain      = "TRB",
    pair_file  = file.path(topdir, "TRBV_TRBJ.csv"),
    cdr3_baseline = cdr3_baseline,
    output_file = file.path(topdir, "TRBV_TRBJ_cdr3.csv")
  )

  # Real alpha + dummy beta
  df_alpha_dummy <- pair_alpha_with_dummy_beta(
    df_A        = df_A,
    output_file = file.path(topdir, "chainA_dummyB_pair.csv")
  )
  prepare_epitope_input(
    df_paired   = df_alpha_dummy,
    peptide     = peptide,
    mhc         = mhc,
    output_file = file.path(topdir, paste0("chainA_dummyB_pair_", peptide, ".csv"))
  )

  # Dummy alpha + real beta
  df_beta_dummy <- pair_beta_with_dummy_alpha(
    df_B        = df_B,
    output_file = file.path(topdir, "dummyA_chainB_pair.csv")
  )
  prepare_epitope_input(
    df_paired   = df_beta_dummy,
    peptide     = peptide,
    mhc         = mhc,
    output_file = file.path(topdir, paste0("dummyA_chainB_pair_", peptide, ".csv"))
  )

  invisible(topdir)
}


### Main ###

input_dir       <- "/Users/roessner/Documents/PostDoc/Data/MixTCRviz/data_raw/HomoSapiens"
base_output_dir <- "../dummy_chains"

baseline      <- MixTCRviz::baseline_HomoSapiens
cdr3_baseline <- baseline$countCDR3.VJL

dico <- list(
  "LLWNGPMAV"  = "HLA_A0201",
  "ELAGIGILTV" = "HLA_A0201",
  "GILGFVFTL"  = "HLA_A0201"
)

for (peptide in names(dico)) {
  mhc <- dico[[peptide]]
  message("Running pipeline for ", peptide, " / ", mhc)
  run_pipeline_for_peptide(
    peptide         = peptide,
    mhc             = mhc,
    base_output_dir = base_output_dir,
    input_dir       = input_dir,
    cdr3_baseline   = cdr3_baseline
  )
}
