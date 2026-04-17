# Computing motifs for AF3 benchmark
# Rita Ann Roessner (ritaann.roessner@unil.ch)
# 2026-04-14

library(MixTCRviz)
library(ggpubr)

make_info <- function(gene, input_name = "Input", baseline_name = "Baseline", model = "Model_default") {
  x <- c(gene, input_name, baseline_name, model)
  names(x) <- c("gene", "input1.name", "baseline.name", "model")
  x
}

baseline.model <- MixTCRviz::baseline_HomoSapiens

df = read.csv('AF3_class_I_motif.csv')
models = unique(df$model)
print(models)
df_factors = read.csv('AF3_class_I_factor.csv')

iptm_threshold = 0.8
out_path = 'motifs/all'

df_filtered <- df[df$AF3_iptm_pair_mean > iptm_threshold , ]
motif <- MixTCRviz(input1=df_filtered, output.path = out_path)

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

str(total_countV[["TRA"]])

ch <- "TRA"
info <- make_info(paste0(ch, "V"))
pVA <- plotVJ(count.es=total_countV[[ch]], count.rep=baseline.model$countV[[ch]], sd.rep = baseline.model$sdV[[ch]], info = info)
info <- make_info(paste0(ch, "J"))
pJA <- plotVJ(count.es=total_countJ[[ch]], count.rep=baseline.model$countJ[[ch]], sd.rep = baseline.model$sdJ[[ch]], info = info)

ch <- "TRB"
info <- make_info(paste0(ch, "V"), model=model)
pVB <- plotVJ(count.es=total_countV[[ch]], count.rep=baseline.model$countV[[ch]], sd.rep = baseline.model$sdV[[ch]], info = info)
info <- make_info(paste0(ch, "J"), model=model)
pJB <- plotVJ(count.es=total_countJ[[ch]], count.rep=baseline.model$countJ[[ch]], sd.rep = baseline.model$sdJ[[ch]], info = info)

p_all <- ggarrange(pVA, pJA, pVB, pJB, ncol = 4, nrow = 1)
ggsave("weighted_motif.png", plot = p_all, width = 16, height = 4, dpi = 300)
