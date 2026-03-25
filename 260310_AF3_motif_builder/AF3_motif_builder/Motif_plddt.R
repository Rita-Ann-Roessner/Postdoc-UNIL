library(MixTCRviz)
library(ggplot2)
library(dplyr)


### Helper functions ###
source("plotting_functions.R")

# helper to parse a cell that contains a list-like string of numbers
parse_num_list <- function(x) {
  if (is.null(x) || is.na(x) || x == "") return(numeric(0))
  
  x <- gsub("\\[|\\]", "", x)
  vals <- unlist(strsplit(x, "[,[:space:]]+"))
  vals <- vals[nzchar(vals)]
  as.numeric(vals)
}

# mean of all plddt values across cdr columns for one row
row_plddt_mean <- function(...) {
  inputs <- list(...)
  
  vals <- unlist(lapply(inputs, parse_num_list))
  
  if (length(vals) == 0) return(NA_real_)
  mean(vals, na.rm = TRUE)
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

infoV <- c("TRAV", "Input", "Baseline", "Model_default")
names(infoV) <- c("gene", "input1.name", "baseline.name", "model")

# load plddt csv
df_plddt = read.csv(file.path(topdir, 'cdr_plddts.csv'))

# add plddt to stats
plddt_vec <- sapply(split(df_plddt, df_plddt$TRAV), function(d) {
  vals <- c(
    unlist(lapply(d$CDR1A, parse_num_list)),
    unlist(lapply(d$CDR2A, parse_num_list)),
    unlist(lapply(d$CDR3AV, parse_num_list))
  )
  mean(vals, na.rm = TRUE)
})

if (is.null(es$plddt)) es$plddt <- list()
es$plddt[["TRA"]] <- plddt_vec

pV <- plotVJ(
  count.es = es$countV[["TRA"]],
  count.rep = baseline.stat$countV[["TRA"]],
  plddt.es = es$plddt[["TRA"]],
  info = infoV,
  species = es$species
)

print(pV)


