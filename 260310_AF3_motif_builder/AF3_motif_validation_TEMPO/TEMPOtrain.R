library(TEMPOtrain)
library(MixTCRviz)

# train on SEQTR data from Liu et al. or on AF3 motif builder top predictions 

batch <- c('SEQTR')
  
predictor <- TEMPOtrain(input.train = paste0("A0201_LLWNGPMAV_", batch,".csv"), output.path = batch,
            input.pred = "validation.csv", filename.train = "A0201_LLWNGPMAV", filename.pred="A0201_LLWNGPMAV_pred", 
            build.prank = T, compute.prank=T, write.data.pred = T, write.data.train = T, write.predictor = T)


batch <- c('AF3')

predictor <- TEMPOtrain(input.train = paste0("A0201_LLWNGPMAV_", batch, ".csv"), output.path  = paste0(batch, "_default_baseline"),
                        input.pred   = "validation.csv", filename.train = "A0201_LLWNGPMAV", filename.pred  = "A0201_LLWNGPMAV_pred",
                        build.prank = T, compute.prank = T, write.data.pred = T, write.data.train = T, write.predictor = T)


# generate baseline
bg <- read.csv(paste0("A0201_LLWNGPMAV_baseline_", batch,".csv"))
baseline <- build_stat(input = bg)  
saveRDS(baseline, paste0("A0201_LLWNGPMAV_baseline_", batch,".rds"))

# prediction with custom baseline
predictor <- TEMPOtrain(input.train = paste0("A0201_LLWNGPMAV_", batch, ".csv"), output.path  = batch,
  input.pred   = "validation.csv", baseline = paste0("A0201_LLWNGPMAV_baseline_", batch,".rds"), filename.train = "A0201_LLWNGPMAV", filename.pred  = "A0201_LLWNGPMAV_pred",
  build.prank = T, compute.prank = T, write.data.pred = T, write.data.train = T, write.predictor = T)


# apparently I can save the predictor (with the built prank and then load it into input.pred=predictor$predictor)