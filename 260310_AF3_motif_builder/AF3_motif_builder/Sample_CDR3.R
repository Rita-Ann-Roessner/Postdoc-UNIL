library(MixTCRviz)

# Function to draw a CDR3 for a given chain and V/J combination
draw_random_cdr3 <- function(chain, v_seg, j_seg, cdr3_baseline) {
  
  # chain: "TRA" or "TRB"
  if (!chain %in% c("TRA", "TRB")) stop("chain must be 'TRA' or 'TRB'")
  
  # Check if V/J exists in baseline
  key <- paste0(v_seg, "_", j_seg)
  if (!key %in% names(cdr3_baseline[[chain]])) {
    return(NA)  # return NA if missing
  }
  
  vj_counts <- cdr3_baseline[[chain]][[key]]
  
  # Generate all length names available for this V/J
  lengths_available <- names(vj_counts)
  
  # Sum counts per length
  lengths_cnt <- sapply(lengths_available, function(len_name) {
    sum(vj_counts[[len_name]][,1], na.rm = TRUE)
  })
  
  # Compute fractions
  fractions <- lengths_cnt / sum(lengths_cnt)
  
  # Draw one random length
  random_length <- sample(lengths_available, size = 1, prob = fractions)
  
  # Get count matrix for that length
  count_mat <- as.matrix(vj_counts[[random_length]])
  
  # Compute position-wise probabilities
  prob_mat <- apply(count_mat, 2, function(col) col / sum(col))
  amino_acids <- rownames(prob_mat)
  L <- ncol(prob_mat)
  
  # Draw one sequence
  seq_vec <- sapply(1:L, function(pos) sample(amino_acids, 1, prob = prob_mat[, pos]))
  paste0(seq_vec, collapse = "")
}

# Load baseline
baseline <- MixTCRviz::baseline_HomoSapiens
cdr3_baseline <- baseline$countCDR3.VJL

# TRA
df <- read.csv('../TRAV_TRAJ.csv')
df$random_CDR3 <- mapply(
  draw_random_cdr3,
  chain = "TRA",
  v_seg = df$TRAV,
  j_seg = df$TRAJ,
  MoreArgs = list(cdr3_baseline = cdr3_baseline)
)
df <- na.omit(df)
write.csv(df, '../TRAV_TRAJ_cdr3.csv', row.names = FALSE)

#TRB
df <- read.csv('../TRBV_TRBJ.csv')
df$random_CDR3 <- mapply(
  draw_random_cdr3,
  chain = "TRB",
  v_seg = df$TRBV,
  j_seg = df$TRBJ,
  MoreArgs = list(cdr3_baseline = cdr3_baseline)
)
df <- na.omit(df)
write.csv(df, '../TRBV_TRBJ_cdr3.csv', row.names = FALSE)