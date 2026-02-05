library(MixTCRviz)

setwd("/Users/roessner/Documents/PostDoc/Data/260202_pca_exp_TCRs")
getwd()

df <- read.csv("data_julien/LAU5013/YF_LAU5013_sc_WT/YF_LAU5013_sc_WT_AF3.tsv", sep="\t")

# For alpha chain
df$cdr1_TRA <- sapply(df$TRAV, function(v) {
  MixTCRviz::cdr123$HomoSapiens$TRA[v, "CDR1"]
})

df$cdr2_TRA <- sapply(df$TRAV, function(v) {
  MixTCRviz::cdr123$HomoSapiens$TRA[v, "CDR2"]
})

# For beta chain
df$cdr1_TRB <- sapply(df$TRBV, function(v) {
  MixTCRviz::cdr123$HomoSapiens$TRB[v, "CDR1"]
})

df$cdr2_TRB <- sapply(df$TRBV, function(v) {
  MixTCRviz::cdr123$HomoSapiens$TRB[v, "CDR2"]
})

write.csv(df, "data_julien/LAU5013/YF_LAU5013_sc_WT/YF_LAU5013_sc_WT_AF3_anno.csv", row.names = FALSE)