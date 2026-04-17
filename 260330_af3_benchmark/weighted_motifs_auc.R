# Computing motifs for AF3 benchmark
# Rita Ann Roessner (ritaann.roessner@unil.ch)
# 2026-04-17

library(MixTCRviz)
library(ggpubr)

make_info <- function(gene, input_name = "Input", baseline_name = "Baseline", model = "Model_default") {
  x <- c(gene, input_name, baseline_name, model)
  names(x) <- c("gene", "input1.name", "baseline.name", "model")
  x
}

df = read.csv('AF3_class_I_motif.csv')
models = unique(df$model)
print(models)
df_factors = read.csv('AF3_class_I_factor.csv')

thresholds <- c(0.6, 1)

for (t in thresholds) {
  df_filtered <- df[df$auc < t , ]
  motif <- MixTCRviz(input1=df_filtered)

  total_countV <- list(TRA = NULL, TRB = NULL)
  total_countJ <- list(TRA = NULL, TRB = NULL)
  
  for (model in models) {
    factor <- df_factors$factor[df_factors$model == model][1]
    es <- motif$stat[[model]]
    
    for (ch in c("TRA", "TRB")) {
      currentV <- es$countV[[ch]] * factor
      currentJ <- es$countJ[[ch]] * factor
      
      if (is.null(total_countV[[ch]])) {
        total_countV[[ch]] <- currentV
      } else {
        all_names <- union(names(total_countV[[ch]]), names(currentV))
        tmp1 <- setNames(rep(0, length(all_names)), all_names)
        tmp2 <- tmp1
        tmp1[names(total_countV[[ch]])] <- total_countV[[ch]]
        tmp2[names(currentV)] <- currentV
        total_countV[[ch]] <- tmp1 + tmp2
      }
      
      if (is.null(total_countJ[[ch]])) {
        total_countJ[[ch]] <- currentJ
      } else {
        all_names <- union(names(total_countJ[[ch]]), names(currentJ))
        tmp1 <- setNames(rep(0, length(all_names)), all_names)
        tmp2 <- tmp1
        tmp1[names(total_countJ[[ch]])] <- total_countJ[[ch]]
        tmp2[names(currentJ)] <- currentJ
        total_countJ[[ch]] <- tmp1 + tmp2
      }
    }
  }
  if (t == 0.6) {
    input1_total_countV <- total_countV
    input1_total_countJ <- total_countJ
  } else {
    input2_total_countV <- total_countV
    input2_total_countJ <- total_countJ
  }
  
}


ch <- "TRA"
info <- make_info(paste0(ch, "V"), input_name = "TCRs (AUC < 0.6)", baseline_name = "Cognate TCRs")
pVA <- plotVJ(count.es=input1_total_countV[[ch]], count.rep=input2_total_countV[[ch]], info = info)
info <- make_info(paste0(ch, "J"), input_name = "TCRs (AUC < 0.6)", baseline_name = "Cognate TCRs")
pJA <- plotVJ(count.es=input1_total_countJ[[ch]], count.rep=input2_total_countJ[[ch]], info = info)

ch <- "TRB"
info <- make_info(paste0(ch, "V"), input_name = "TCRs (AUC < 0.6)", baseline_name = "Cognate TCRs")
pVB <- plotVJ(count.es=input1_total_countV[[ch]], count.rep=input2_total_countV[[ch]], info = info)
info <- make_info(paste0(ch, "J"), input_name = "TCRs (AUC < 0.6)", baseline_name = "Cognate TCRs")
pJB <- plotVJ(count.es=input1_total_countJ[[ch]], count.rep=input2_total_countJ[[ch]], info = info)

p_all <- ggarrange(pVA, pJA, pVB, pJB, ncol = 4, nrow = 1)
print(p_all)
ggsave(paste0("motif_low_auc_class_I.png"), plot = p_all, width = 16, height = 4, dpi = 300)



