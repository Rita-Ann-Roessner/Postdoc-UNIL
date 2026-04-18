library(MixTCRviz)
library(dplyr)
library(tidyr)

### Helper Functions ###

compute_enrichment <- function(baseline, models, n = NULL, cols = c("TRAV", "TRAJ", "TRBV", "TRBJ")) {
  count_df <- function(df, source_label) {
    df %>%
      select(all_of(cols)) %>%
      pivot_longer(everything(), names_to = "column", values_to = "value") %>%
      group_by(column, value) %>%
      summarise(count = n(), .groups = "drop") %>%
      rename(!!paste0(source_label, "_count") := count)
  }
  
  baseline_counts <- count_df(baseline, "baseline")
  models_counts   <- count_df(models, "models")
  
  final_df <- baseline_counts %>%
    left_join(models_counts, by = c("column", "value")) %>%
    mutate(models_count = replace_na(models_count, 0),
           ratio = models_count / baseline_count)
  
  if (!is.null(n)) {
    final_df <- final_df %>%
      group_by(column) %>%
      slice_max(order_by = ratio, n = n, with_ties = FALSE) %>%
      ungroup()
  }
  
  return(final_df)
}


draw_random_cdr3_multi <- function(chain, v_seg, j_seg, cdr3_baseline, n = 5, len_range = 10:20) {
  if (!chain %in% c("TRA", "TRB")) stop("chain must be 'TRA' or 'TRB'")
  
  key <- paste0(v_seg, "_", j_seg)
  if (!key %in% names(cdr3_baseline[[chain]])) {
    return(NA)
  }
  
  vj_counts <- cdr3_baseline[[chain]][[key]]
  
  lengths_available <- names(vj_counts)
  numeric_lengths <- as.numeric(sub("L_", "", lengths_available))
  valid_lengths <- intersect(len_range, numeric_lengths)
  if (length(valid_lengths) == 0) return(NA)
  
  valid_lengths <- paste0("L_", valid_lengths)
  selected_lengths <- sample(valid_lengths, size = min(n, length(valid_lengths)))
  
  sequences <- c()
  
  for (L in selected_lengths) {
    count_mat <- as.matrix(vj_counts[[L]])
    prob_mat <- apply(count_mat, 2, function(col) col / sum(col))
    amino_acids <- rownames(prob_mat)
    
    seq_vec <- sapply(1:ncol(prob_mat), function(pos) {
      sample(amino_acids, 1, prob = prob_mat[, pos])
    })
    
    sequences <- c(sequences, paste0(seq_vec, collapse = ""))
  }
  
  return(sequences)
}


load_enrichment <- function(peptide, base_dir, threshold = 0.5, n = 10) {
  step1_dir <- file.path(base_dir, paste0("step1_", peptide))
  
  df <- read.table(
    file.path(step1_dir, paste0("chainA_B_random_pair_", peptide, "_output.txt")),
    header = TRUE
  )[, c("id", "AF3_iptm_pair_mean")]
  
  anno <- read.table(
    file.path(step1_dir, paste0("chainA_B_random_pair_", peptide, "_input.txt")),
    header = TRUE
  )
  
  df <- merge(df, anno, by = "id")
  
  models <- df[df$AF3_iptm_pair_mean > threshold, ]
  
  compute_enrichment(df, models, n = n)
}


generate_enriched_vj_pairs <- function(enriched_df, output_dir) {
  for (chain in c("A", "B")) {
    genes <- c(paste0("TR", chain, "V"), paste0("TR", chain, "J"))
    
    lst <- list(
      enriched_df$value[enriched_df$column == genes[1]],
      enriched_df$value[enriched_df$column == genes[2]]
    )
    
    pairs <- expand.grid(lst[[1]], lst[[2]])
    colnames(pairs) <- genes
    
    write.csv(pairs, file.path(output_dir, paste0(genes[1], "_", genes[2], ".csv")), row.names = FALSE)
  }
}


sample_chain_cdr3_multi <- function(chain, pair_file, cdr3_baseline, output_file, n = 5) {
  df <- read.csv(pair_file)
  
  if (chain == "TRA") {
    v_col <- "TRAV"; j_col <- "TRAJ"
  } else if (chain == "TRB") {
    v_col <- "TRBV"; j_col <- "TRBJ"
  } else {
    stop("chain must be 'TRA' or 'TRB'")
  }
  
  results <- mapply(
    draw_random_cdr3_multi,
    chain = chain,
    v_seg = df[[v_col]],
    j_seg = df[[j_col]],
    MoreArgs = list(cdr3_baseline = cdr3_baseline, n = n),
    SIMPLIFY = FALSE
  )
  
  df_expanded <- data.frame(
    V = rep(df[[v_col]], lengths(results)),
    J = rep(df[[j_col]], lengths(results)),
    random_CDR3 = unlist(results)
  )
  colnames(df_expanded)[1:2] <- c(v_col, j_col)
  
  df_expanded <- na.omit(df_expanded)
  write.csv(df_expanded, output_file, row.names = FALSE)
  
  df_expanded
}


pair_alpha_beta_multi <- function(file_a, file_b, output_file, n = 5) {
  df_A <- read.csv(file_a)
  df_A <- df_A[rep(1:nrow(df_A), times = n), ]
  colnames(df_A)[colnames(df_A) == "random_CDR3"] <- "cdr3_TRA"
  
  df_B <- read.csv(file_b)
  df_B <- df_B[rep(1:nrow(df_B), times = n), ]
  colnames(df_B)[colnames(df_B) == "random_CDR3"] <- "cdr3_TRB"
  
  n_A <- nrow(df_A)
  n_B <- nrow(df_B)
  n_max <- max(n_A, n_B)
  
  df_A_expanded <- df_A[sample(1:n_A, n_max, replace = (n_A < n_max)), ]
  df_B_expanded <- df_B[sample(1:n_B, n_max, replace = (n_B < n_max)), ]
  
  df_A_expanded <- df_A_expanded[sample(1:n_max), ]
  df_B_expanded <- df_B_expanded[sample(1:n_max), ]
  
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


run_pipeline_for_peptide <- function(peptide, mhc, base_output_dir, cdr3_baseline,
                                     threshold = 0.5, n_enrichment = 20, n_cdr3 = 5, n_pairs = 2) {
  enriched <- load_enrichment(peptide, base_output_dir, threshold = threshold, n = n_enrichment)

  topdir <- file.path(base_output_dir, paste0("step2_", peptide))
  dir.create(topdir, showWarnings = FALSE, recursive = TRUE)
  
  generate_enriched_vj_pairs(enriched, topdir)
  
  sample_chain_cdr3_multi(
    chain = "TRA",
    pair_file = file.path(topdir, "TRAV_TRAJ.csv"),
    cdr3_baseline = cdr3_baseline,
    output_file = file.path(topdir, "TRAV_TRAJ_cdr3.csv"),
    n = n_cdr3
  )
  
  sample_chain_cdr3_multi(
    chain = "TRB",
    pair_file = file.path(topdir, "TRBV_TRBJ.csv"),
    cdr3_baseline = cdr3_baseline,
    output_file = file.path(topdir, "TRBV_TRBJ_cdr3.csv"),
    n = n_cdr3
  )
  
  df_paired <- pair_alpha_beta_multi(
    file_a = file.path(topdir, "TRAV_TRAJ_cdr3.csv"),
    file_b = file.path(topdir, "TRBV_TRBJ_cdr3.csv"),
    output_file = file.path(topdir, "chainA_B_random_pair.csv"),
    n = n_pairs
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
    cdr3_baseline = cdr3_baseline
  )
}