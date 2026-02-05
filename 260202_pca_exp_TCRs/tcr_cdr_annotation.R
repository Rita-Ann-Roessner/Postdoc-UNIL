library(MixTCRviz)

setwd("/Users/roessner/Documents/PostDoc/Data/260202_pca_exp_TCRs")

df <- read.csv("data_julien/LAU5013/YF_LAU5013_sc_WT/YF_LAU5013_sc_WT_AF3.tsv", sep="\t")

# Filter valid V genes from MixTCRviz
pA <- which(!grepl("*", rownames(MixTCRviz::cdr123$HomoSapiens$TRA), fixed=TRUE))
pB <- which(!grepl("*", rownames(MixTCRviz::cdr123$HomoSapiens$TRB), fixed=TRUE) & 
              rownames(MixTCRviz::cdr123$HomoSapiens$TRB) != "TRBV6-3" &
              rownames(MixTCRviz::cdr123$HomoSapiens$TRB) != "TRBV6-2/6-3" &
              rownames(MixTCRviz::cdr123$HomoSapiens$TRB) != "TRBV12-3/12-4")

# Build lookup tables (named vectors)
cdr1 <- list()
cdr1[["TRA"]] <- setNames(gsub("g","", MixTCRviz::cdr123$HomoSapiens$TRA$CDR1[pA]), rownames(MixTCRviz::cdr123$HomoSapiens$TRA)[pA])
cdr1[["TRB"]] <- setNames(gsub("g","", MixTCRviz::cdr123$HomoSapiens$TRB$CDR1[pB]), rownames(MixTCRviz::cdr123$HomoSapiens$TRB)[pB])

cdr2 <- list()
cdr2[["TRA"]] <- setNames(gsub("g","", MixTCRviz::cdr123$HomoSapiens$TRA$CDR2[pA]), rownames(MixTCRviz::cdr123$HomoSapiens$TRA)[pA])
cdr2[["TRB"]] <- setNames(gsub("g","", MixTCRviz::cdr123$HomoSapiens$TRB$CDR2[pB]), rownames(MixTCRviz::cdr123$HomoSapiens$TRB)[pB])

# Annotate dataframe using the cleaned lookups
df$cdr1_TRA <- sapply(df$TRAV, function(v) if(v %in% names(cdr1$TRA)) cdr1$TRA[v] else NA)
df$cdr2_TRA <- sapply(df$TRAV, function(v) if(v %in% names(cdr2$TRA)) cdr2$TRA[v] else NA)

df$cdr1_TRB <- sapply(df$TRBV, function(v) if(v %in% names(cdr1$TRB)) cdr1$TRB[v] else NA)
df$cdr2_TRB <- sapply(df$TRBV, function(v) if(v %in% names(cdr2$TRB)) cdr2$TRB[v] else NA)

# Save the annotated dataframe
write.csv(df, "data_julien/LAU5013/YF_LAU5013_sc_WT/YF_LAU5013_sc_WT_AF3_anno.csv", row.names = FALSE)
