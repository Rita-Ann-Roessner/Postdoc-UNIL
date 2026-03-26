library(MixTCRviz)
library(ggplot2)
library(dplyr)


### Helper functions ###
source("plotting_functions.R")

make_info <- function(gene, input_name = "Input", baseline_name = "Baseline", model = "Model_default") {
  x <- c(gene, input_name, baseline_name, model)
  names(x) <- c("gene", "input1.name", "baseline.name", "model")
  x
}

# helper to parse a cell that contains a list-like string of numbers
parse_num_list <- function(x) {
  if (is.null(x) || is.na(x) || x == "") return(numeric(0))
  
  x <- gsub("\\[|\\]", "", x)
  vals <- unlist(strsplit(x, "[,[:space:]]+"))
  vals <- vals[nzchar(vals)]
  as.numeric(vals)
}

# compute per gene plddt
make_plddt <- function(gene_col, cols) {
  sapply(split(df_plddt, df_plddt[[gene_col]]), function(d) {
    mean(unlist(lapply(cols, function(col) unlist(lapply(d[[col]], parse_num_list)))), na.rm = TRUE)
  })
}

# Amino-acid alphabet in the same order used by MixTCRviz
aa.list <- c("A","C","D","E","F","G","H","I","K","L",
             "M","N","P","Q","R","S","T","V","W","Y")
N.aa <- length(aa.list)

# Build per-length amino-acid/position pLDDT matrices for one chain
# seq_col: column with CDR3 sequence, e.g. "CDR3A_seq"
# plddt_col: column with per-residue pLDDT list-string, e.g. "CDR3A"
make_plddt_cdr3_matrices <- function(df, seq_col, plddt_col, aa_order = aa.list) {
  
  # keep only rows with usable sequence + plddt
  keep <- !is.na(df[[seq_col]]) & df[[seq_col]] != "" &
    !is.na(df[[plddt_col]]) & df[[plddt_col]] != ""
  df <- df[keep, , drop = FALSE]
  
  if (nrow(df) == 0) return(list())
  
  # parse inputs
  seqs <- as.character(df[[seq_col]])
  plddts <- lapply(df[[plddt_col]], parse_num_list)
  lens <- nchar(seqs)
  
  # keep only rows where sequence length matches number of plddt values
  ok <- mapply(function(s, p) nchar(s) == length(p), seqs, plddts)
  seqs <- seqs[ok]
  plddts <- plddts[ok]
  lens <- lens[ok]
  
  if (length(seqs) == 0) return(list())
  
  out <- list()
  
  for (L in sort(unique(lens))) {
    idx <- which(lens == L)
    seqs_L <- seqs[idx]
    plddts_L <- plddts[idx]
    
    # sum of pLDDT values for each AA x position
    sum_mat <- matrix(0, nrow = length(aa_order), ncol = L,
                      dimnames = list(aa_order, paste0("pos", seq_len(L))))
    
    # number of observations contributing to each AA x position
    n_mat <- matrix(0, nrow = length(aa_order), ncol = L,
                    dimnames = list(aa_order, paste0("pos", seq_len(L))))
    
    for (i in seq_along(seqs_L)) {
      aa_vec <- strsplit(seqs_L[i], "", fixed = TRUE)[[1]]
      p_vec <- plddts_L[[i]]
      
      for (pos in seq_len(L)) {
        aa <- aa_vec[pos]
        if (!is.na(aa) && aa %in% aa_order && !is.na(p_vec[pos])) {
          sum_mat[aa, pos] <- sum_mat[aa, pos] + p_vec[pos]
          n_mat[aa, pos] <- n_mat[aa, pos] + 1
        }
      }
    }
    
    avg_mat <- sum_mat
    avg_mat[n_mat > 0] <- sum_mat[n_mat > 0] / n_mat[n_mat > 0]
    avg_mat[n_mat == 0] <- 0
    
    out[[paste0("L_", L)]] <- avg_mat
  }
  
  out
}


# compute average pLDDT per CDR3 length from a column like CDR3A or CDR3B
make_plddt_by_length <- function(col, prefix = "L_") {
  vals_list <- lapply(df_plddt[[col]], parse_num_list)
  
  # length of each CDR3
  cdr3_len <- sapply(vals_list, length)
  
  # mean pLDDT of each CDR3
  cdr3_mean <- sapply(vals_list, function(x) {
    if (length(x) == 0) return(NA_real_)
    mean(x, na.rm = TRUE)
  })
  
  # keep only non-empty entries
  keep <- cdr3_len > 0 & !is.na(cdr3_mean)
  cdr3_len <- cdr3_len[keep]
  cdr3_mean <- cdr3_mean[keep]
  
  # average pLDDT by length
  out <- tapply(cdr3_mean, cdr3_len, mean, na.rm = TRUE)
  
  out <- as.numeric(out)
  names(out) <- paste0(prefix, names(tapply(cdr3_mean, cdr3_len, mean, na.rm = TRUE)))
  
  out
}

### Annotate motif with AF3 pLDDT values ### 

# 1. Compute stats
topdir <- "../step2_LLWNGPMAV"

df <- read.table(file.path(topdir, "chainA_B_random_pair_LLWNGPMAV_output.txt"), header = TRUE)[, c("id", "AF3_iptm_pair_mean")]
anno <- read.table(file.path(topdir, "chainA_B_random_pair_LLWNGPMAV_input.txt"), header = TRUE)

df <- merge(df, anno, by = "id")

baseline <- df

threshold <- 0.8
model <- df[df$AF3_iptm_pair_mean > threshold, ]

motif <- MixTCRviz(input1=model, input2=baseline) #plot = FALSE
es <- motif$stat$Model_default

# 2. Plot motif
es <- motif$stat[[1]]
baseline.stat <- MixTCRviz::build_stat(baseline, chain = "AB", species = es$species, comp.VJL = 0)

# load plddt csv
df_plddt = read.csv(file.path(topdir, 'cdr_plddts.csv'))

# add plddt to VJ stats
if (is.null(es$plddt)) es$plddt <- list()
es$plddtV[["TRA"]] <- make_plddt("TRAV", c("CDR1A", "CDR2A", "CDR3AV"))
es$plddtV[["TRB"]] <- make_plddt("TRBV", c("CDR1B", "CDR2B", "CDR3BV"))
es$plddtJ[["TRA"]] <- make_plddt("TRAJ", c("CDR3AJ"))
es$plddtJ[["TRB"]] <- make_plddt("TRBJ", c("CDR3BJ"))

pVA <- plotVJ(
  count.es = es$countV[["TRA"]],
  count.rep = baseline.stat$countV[["TRA"]],
  plddt.es = es$plddtV[["TRA"]],
  info = make_info("TRAV"),
  species = es$species,
  show.plddt.legend = TRUE
)

pVB <- plotVJ(
  count.es = es$countV[["TRB"]],
  count.rep = baseline.stat$countV[["TRB"]],
  plddt.es = es$plddtV[["TRB"]],
  info = make_info("TRBV"),
  species = es$species
)

pJA <- plotVJ(
  count.es = es$countJ[["TRA"]],
  count.rep = baseline.stat$countJ[["TRA"]],
  plddt.es = es$plddtJ[["TRA"]],
  info = make_info("TRAJ"),
  species = es$species
)

pJB <- plotVJ(
  count.es = es$countJ[["TRB"]],
  count.rep = baseline.stat$countJ[["TRB"]],
  plddt.es = es$plddtJ[["TRB"]],
  info = make_info("TRBJ"),
  species = es$species
)

# 3. Build the remaining panels for the final motif figure
#    (length distribution + CDR3 logos)

infoA <- c(
  chain = "A",
  input1.name = "Input",
  baseline.name = "Baseline",
  model = "Model_default"
)

infoB <- c(
  chain = "B",
  input1.name = "Input",
  baseline.name = "Baseline",
  model = "Model_default"
)

# add plddt to length distribution plots
es$plddtL[["TRA"]] <- make_plddt_by_length("CDR3A")
es$plddtL[["TRB"]]  <- make_plddt_by_length("CDR3B")

# Length-distribution plots
pLDA <- plotLD(
  countL.es  = es$countL[["TRA"]],
  countL.rep = baseline.stat$countL[["TRA"]],
  plddtL.es = es$plddtL[["TRA"]],
  info = infoA
)

pLDB <- plotLD(
  countL.es  = es$countL[["TRB"]],
  countL.rep = baseline.stat$countL[["TRB"]],
  plddtL.es = es$plddtL[["TRB"]],
  info = infoB
)

# add plddt to CDR3 logo objects
es$plddtCDR3.L[["TRA"]] <- make_plddt_cdr3_matrices(
  df = df_plddt,
  seq_col = "CDR3A_seq",
  plddt_col = "CDR3A"
)

es$plddtCDR3.L[["TRB"]] <- make_plddt_cdr3_matrices(
  df = df_plddt,
  seq_col = "CDR3B_seq",
  plddt_col = "CDR3B"
)

# CDR3 logo objects
CDR3A <- plotCDR3(
  countL.es     = es$countL[["TRA"]],
  countL.rep    = baseline.stat$countL[["TRA"]],
  countCDR3.es  = es$countCDR3.L[["TRA"]],
  countCDR3.rep = baseline.stat$countCDR3.L[["TRA"]],
  plddtCDR3.L.es = es$plddtCDR3.L[["TRA"]],
  info = infoA
)

CDR3B <- plotCDR3(
  countL.es     = es$countL[["TRB"]],
  countL.rep    = baseline.stat$countL[["TRB"]],
  countCDR3.es  = es$countCDR3.L[["TRB"]],
  countCDR3.rep = baseline.stat$countCDR3.L[["TRB"]],
  plddtCDR3.L.es = es$plddtCDR3.L[["TRB"]],
  info = infoB
)

# combine input and baseline CDR3 motifs for each chain

cdr3_panel_A <- ggpubr::ggarrange(
  CDR3A$ES_max,
  CDR3A$Baseline_max,
  nrow = 2
)

cdr3_panel_B <- ggpubr::ggarrange(
  CDR3B$ES_max,
  CDR3B$Baseline_max,
  nrow = 2
)

pg.all <- list()

pg.all[["TRA"]] <- ggpubr::ggarrange(
  pVA, pJA,
  pLDA, cdr3_panel_A,
  ncol = 2, nrow = 2
)

pg.all[["TRB"]] <- ggpubr::ggarrange(
  pVB, pJB,
  pLDB, cdr3_panel_B,
  ncol = 2, nrow = 2
)

fig <- ggpubr::ggarrange(
  pg.all[["TRA"]],
  pg.all[["TRB"]],
  ncol = 2
)

grid::grid.newpage()
grid::grid.draw(fig)

ggplot2::ggsave(
  "final_motif.png",
  fig,
  width = 20,
  height = 10,
  dpi = 300
)

