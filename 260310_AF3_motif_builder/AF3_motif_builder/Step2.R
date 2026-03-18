library(MixTCRviz)
library(dplyr)
library(tidyr)

### Helper functions ###
compute_enrichment <- function(baseline, models, n = NULL, cols = c("TRAV", "TRAJ", "TRBV", "TRBJ")) {
  # Function to count and label source
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
  
  # Join counts and compute ratio
  final_df <- baseline_counts %>%
    left_join(models_counts, by = c("column", "value")) %>%
    mutate(models_count = replace_na(models_count, 0),
           ratio = models_count / baseline_count)
  
  # select top n
  if (!is.null(n)) {
    final_df <- final_df %>%
      group_by(column) %>%                # remove this line if you want global top n instead
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
  
  # Only lengths available in this V/J
  lengths_available <- names(vj_counts)
  numeric_lengths <- as.numeric(sub("L_", "", lengths_available))
  valid_lengths <- intersect(len_range, numeric_lengths)
  if (length(valid_lengths) == 0) return(NA)
  
  valid_lengths <- paste0("L_", valid_lengths)
  
  # Sample up to n lengths
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

### Step 2 ### 
# 1. Compute enrichment for high confidence models
topdir <- "../step1_LLWNGPMAV"

df <- read.table(file.path(topdir, "chainA_B_random_pair_LLWNGPMAV_output.txt"), header = TRUE)[, c("id", "AF3_iptm_pair_mean")]
anno <- read.table(file.path(topdir, "chainA_B_random_pair_LLWNGPMAV_input.txt"), header = TRUE)

df <- merge(df, anno, by = "id")

baseline <- df

threshold <- 0.5
models <- df[df$AF3_iptm_pair_mean > threshold, ]

df = compute_enrichment(baseline, models, n=10)

# 2. Generate V/J pairings for A/B chain
topdir <- "../step2_LLWNGPMAV"
dir.create(topdir, showWarnings = FALSE, recursive = TRUE)

chains <- c("A", "B")

for (chain in chains) {
  genes <- c(paste0("TR", chain, "V"), paste0("TR", chain, "J"))
  lst <- list()
  
  for (i in seq_along(genes)) {
    gene <- genes[i]
    lst[[i]] <- df$value[df$column == gene]
  }

  pairs <- expand.grid(lst[[1]], lst[[2]])
  colnames(pairs) <- genes
  
  write.csv(pairs, file.path(topdir, paste0(genes[1], "_", genes[2], ".csv")), row.names = FALSE)
}

# 3. Sample CDR3s from baseline of different length
baseline <- MixTCRviz::baseline_HomoSapiens
cdr3_baseline <- baseline$countCDR3.VJL

# TRA
df <- read.csv(file.path(topdir, 'TRAV_TRAJ.csv'))
results <- mapply(
  draw_random_cdr3_multi,
  chain = "TRA",
  v_seg = df$TRAV,
  j_seg = df$TRAJ,
  MoreArgs = list(cdr3_baseline = cdr3_baseline, n = 5),
  SIMPLIFY = FALSE
)

df_expanded <- data.frame(
  TRAV = rep(df$TRAV, lengths(results)),
  TRAJ = rep(df$TRAJ, lengths(results)),
  random_CDR3 = unlist(results)
)

df_expanded <- na.omit(df_expanded)
write.csv(df_expanded, file.path(topdir, 'TRAV_TRAJ_cdr3.csv'), row.names = FALSE)

# TRB
df <- read.csv(file.path(topdir, 'TRBV_TRBJ.csv'))
results <- mapply(
  draw_random_cdr3_multi,
  chain = "TRB",
  v_seg = df$TRBV,
  j_seg = df$TRBJ,
  MoreArgs = list(cdr3_baseline = cdr3_baseline, n = 5),
  SIMPLIFY = FALSE
)

df_expanded <- data.frame(
  TRBV = rep(df$TRBV, lengths(results)),
  TRBJ = rep(df$TRBJ, lengths(results)),
  random_CDR3 = unlist(results)
)

df_expanded <- na.omit(df_expanded)
write.csv(df_expanded, file.path(topdir, 'TRBV_TRBJ_cdr3.csv'), row.names = FALSE)

