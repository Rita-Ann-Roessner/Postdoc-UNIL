library(MixTCRviz)

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


pair_alpha_beta <- function(file_a, file_b, output_file, seed = 42) {
  df_A <- read.csv(file_a)
  df_B <- read.csv(file_b)
  
  colnames(df_A)[colnames(df_A) == "random_CDR3"] <- "cdr3_TRA"
  colnames(df_B)[colnames(df_B) == "random_CDR3"] <- "cdr3_TRB"
  
  n_A <- nrow(df_A)
  n_B <- nrow(df_B)
  
  if (n_A == 0 || n_B == 0) {
    stop("One of the chain dataframes is empty after CDR3 sampling.")
  }
  
  n_max <- max(n_A, n_B)
  
  set.seed(seed)
  
  df_A_expanded <- df_A[sample(seq_len(n_A), n_max, replace = (n_A < n_max)), , drop = FALSE]
  df_B_expanded <- df_B[sample(seq_len(n_B), n_max, replace = (n_B < n_max)), , drop = FALSE]
  
  df_A_expanded <- df_A_expanded[sample(seq_len(n_max)), , drop = FALSE]
  df_B_expanded <- df_B_expanded[sample(seq_len(n_max)), , drop = FALSE]
  
  rownames(df_A_expanded) <- NULL
  rownames(df_B_expanded) <- NULL
  
  df_paired <- cbind(df_A_expanded, df_B_expanded)
  write.csv(df_paired, output_file, row.names = FALSE)
  
  df_paired
}


prepare_epitope_input <- function(df_paired, peptide, mhc, species = "HomoSapiens", output_file) {
  df <- df_paired
  df$peptide <- peptide
  df$MHC <- mhc
  df$species <- species
  
  write.csv(df, output_file, row.names = FALSE)
  df
}


run_pipeline_for_peptide <- function(peptide, mhc, base_output_dir, input_dir, cdr3_baseline) {
  topdir <- file.path(base_output_dir, paste0("step1_", peptide))
  dir.create(topdir, showWarnings = FALSE, recursive = TRUE)
  
  generate_vj_pairs("A", input_dir, topdir)
  generate_vj_pairs("B", input_dir, topdir)
  
  sample_chain_cdr3(
    chain = "TRA",
    pair_file = file.path(topdir, "TRAV_TRAJ.csv"),
    cdr3_baseline = cdr3_baseline,
    output_file = file.path(topdir, "TRAV_TRAJ_cdr3.csv")
  )
  
  sample_chain_cdr3(
    chain = "TRB",
    pair_file = file.path(topdir, "TRBV_TRBJ.csv"),
    cdr3_baseline = cdr3_baseline,
    output_file = file.path(topdir, "TRBV_TRBJ_cdr3.csv")
  )
  
  df_paired <- pair_alpha_beta(
    file_a = file.path(topdir, "TRAV_TRAJ_cdr3.csv"),
    file_b = file.path(topdir, "TRBV_TRBJ_cdr3.csv"),
    output_file = file.path(topdir, "chainA_B_random_pair.csv")
  )
  
  prepare_epitope_input(
    df_paired = df_paired,
    peptide = peptide,
    mhc = mhc,
    output_file = file.path(topdir, paste0("chainA_B_random_pair_", peptide, ".csv"))
  )
  
  invisible(topdir)
}


### Main ###

input_dir <- "/Users/roessner/Documents/PostDoc/Data/MixTCRviz/data_raw/HomoSapiens"
base_output_dir <- ".."

baseline <- MixTCRviz::baseline_HomoSapiens
cdr3_baseline <- baseline$countCDR3.VJL

dico <- list(
  "LLWNGPMAV" = "HLA_A0201",
  "ELAGIGILTV" = "HLA_A0201",
  "GILGFVFTL" = "HLA_A0201"
)

for (peptide in names(dico)) {
  mhc <- dico[[peptide]]
  message("Running pipeline for ", peptide, " / ", mhc)
  run_pipeline_for_peptide(
    peptide = peptide,
    mhc = mhc,
    base_output_dir = base_output_dir,
    input_dir = input_dir,
    cdr3_baseline = cdr3_baseline
  )
}