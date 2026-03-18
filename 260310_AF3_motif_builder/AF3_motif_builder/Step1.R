library(MixTCRviz)

### Helper Functions ###
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

### STEP 1 ###
topdir <- "../step1_LLWNGPMAV"
dir.create(topdir, showWarnings = FALSE, recursive = TRUE)

# 1. Generate V/J pairings for A/B chain
chains <- c("A", "B")

for (chain in chains) {
  genes <- c(paste0("TR", chain, "V"), paste0("TR", chain, "J"))
  lst <- list()
  
  for (i in seq_along(genes)) {
    gene <- genes[i]
    df_gene <- read.csv(paste0("/Users/roessner/Documents/PostDoc/Data/MixTCRviz/data_raw/HomoSapiens/", gene, ".csv"))
    
    lst[[i]] <- df_gene[[1]]  # first column
  }
  
  pairs <- expand.grid(lst[[1]], lst[[2]])
  colnames(pairs) <- genes
  
  write.csv(pairs, paste0(genes[1], "_", genes[2], ".csv"), row.names = FALSE)
}

# 2. Sample CDR3s from baseline
# Load baseline
baseline <- MixTCRviz::baseline_HomoSapiens
cdr3_baseline <- baseline$countCDR3.VJL

# TRA
df <- read.csv(file.path(topdir, 'TRAV_TRAJ.csv'))
df$random_CDR3 <- mapply(
  draw_random_cdr3,
  chain = "TRA",
  v_seg = df$TRAV,
  j_seg = df$TRAJ,
  MoreArgs = list(cdr3_baseline = cdr3_baseline)
)
df <- na.omit(df)
write.csv(df, file.path(topdir, 'TRAV_TRAJ_cdr3.csv'), row.names = FALSE)

#TRB
df <- read.csv(file.path(topdir, 'TRBV_TRBJ.csv'))
df$random_CDR3 <- mapply(
  draw_random_cdr3,
  chain = "TRB",
  v_seg = df$TRBV,
  j_seg = df$TRBJ,
  MoreArgs = list(cdr3_baseline = cdr3_baseline)
)
df <- na.omit(df)
write.csv(df, file.path(topdir, 'TRBV_TRBJ_cdr3.csv'), row.names = FALSE)

# 3. Randomly pair A/B chain
df_A <- read.csv(file.path(topdir, "TRAV_TRAJ_cdr3.csv"))
colnames(df_A)[colnames(df_A) == "random_CDR3"] <- "cdr3_TRA"

df_B <- read.csv(file.path(topdir, "TRBV_TRBJ_cdr3.csv"))
colnames(df_B)[colnames(df_B) == "random_CDR3"] <- "cdr3_TRB"

n_A <- nrow(df_A)
n_B <- nrow(df_B)

n_max <- max(n_A, n_B)

set.seed(42)

# Randomly sample rows to match the longer dataframe
df_A_expanded <- df_A[sample(1:n_A, n_max, replace = (n_A < n_max)), ]
df_B_expanded <- df_B[sample(1:n_B, n_max, replace = (n_B < n_max)), ]

# Shuffle both dataframes independently
df_A_expanded <- df_A_expanded[sample(1:n_max), ]
df_B_expanded <- df_B_expanded[sample(1:n_max), ]

rownames(df_A_expanded) <- NULL
rownames(df_B_expanded) <- NULL

df_paired <- cbind(df_A_expanded, df_B_expanded)

write.csv(df_paired, file.path(topdir, "chainA_B_random_pair.csv"), row.names = FALSE)

# 4. Prepare AF3 input for epitopes with known motifs
dico <- list(
  "LLWNGPMAV" = "HLA_A0201",
  "ELAGIGILTV" = "HLA_A0201",
  "GILGFVFTL" = "HLA_A0201"
)

for (peptide in names(dico)) {
  MHC <- dico[[peptide]]
  
  df <- read.csv(file.path(topdir, "chainA_B_random_pair.csv"))
  
  df$peptide <- peptide
  df$MHC <- MHC
  df$species <- "HomoSapiens"
  
  write.csv(df, file.path(topdir, paste0("chainA_B_random_pair_", peptide, ".csv")), row.names = FALSE)
}