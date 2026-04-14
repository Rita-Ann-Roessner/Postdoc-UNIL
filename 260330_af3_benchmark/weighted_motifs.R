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

df = read.csv('test.csv')
models = unique(df$model)
df_factors = read.csv('AF3_class_I_factor.csv')

iptm_threshold = 0.5
out_path = 'motifs/medium_conf'

df_filtered <- df[df$AF3_iptm_pair_mean > iptm_threshold , ]
motif <- MixTCRviz(input1=df_filtered, output.path = out_path)

for (model in models){
  print(model)
  factor = df_factors[df_factors$model == model, "factor"]
  str(factor)
  es <- motif$stat[[model]]
  
  ch <- "TRA"
  es$countV[[ch]] <- es$countV[[ch]] * factor
  es$countJ[[ch]] <- es$countJ[[ch]] * factor
  motif$stat[[model]] <- es
  View(motif)
}


for (model in models){
  print(model)
  factor = df_factors[df_factors$model == model, "factor"]
  es <- motif$stat[[model]]
  
  ch <- "TRA"
  info <- make_info(paste0(ch, "V"), model=model)
  pVA <- plotVJ(count.es=es$countV[[ch]], count.rep=baseline.model$countV[[ch]], info = info)
  info <- make_info(paste0(ch, "J"), model=model)
  pJA <- plotVJ(count.es=es$countJ[[ch]], count.rep=baseline.model$countJ[[ch]], info = info)

  ch <- "TRA"
  info <- make_info(paste0(ch, "V"), model=model)
  pVB <- plotVJ(count.es=es$countV[[ch]], count.rep=baseline.model$countV[[ch]], info = info)
  info <- make_info(paste0(ch, "J"), model=model)
  pJB <- plotVJ(count.es=es$countJ[[ch]], count.rep=baseline.model$countJ[[ch]], info = info)
  
  p_all <- ggarrange(pVA, pJA, pVB, pJB, ncol = 4, nrow = 1)
  print(p_all)
}

