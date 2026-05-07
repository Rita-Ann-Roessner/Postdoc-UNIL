library(TEMPOtrain)
library(MixTCRviz)

# train on SEQTR data from Liu et al. 
batch <- c('SEQTR')
  
predictor <- TEMPOtrain(input.train = paste0("A0201_LLWNGPMAV_", batch,".csv"), output.path = batch,
            input.pred = "validation.csv", filename.train = "A0201_LLWNGPMAV", filename.pred="A0201_LLWNGPMAV_pred", 
            build.prank = T, compute.prank=T, write.data.pred = T, write.data.train = T, write.predictor = T)

# train on top 1% AF3 predictions of vaccinated donor

batch <- c('LAU5013_scPreTCR_AF3')

predictor <- TEMPOtrain(input.train = paste0("A0201_LLWNGPMAV_", batch,".csv"), output.path = batch,
                        input.pred = "validation.csv", filename.train = "A0201_LLWNGPMAV", filename.pred="A0201_LLWNGPMAV_pred", 
                        build.prank = T, compute.prank=T, write.data.pred = T, write.data.train = T, write.predictor = T)

# train on AF3 motif builder top predictions 

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

# dummy chain vs random paring
# dummy chains
peptides <- c('ELAGIGILTV', 'GILGFVFTL', 'LLWNGPMAV')

for (peptide in peptides) {
  for (chain in chains) {
    indir  <- file.path('../dummy_chains', paste0('step1_dummy_', peptide))
    outdir <- file.path(indir, 'TEMPO')
    
    predictor <- TEMPOtrain(input.train = file.path(indir, 'model.csv'), output.path = outdir,
      input.pred      = 'validation.csv', filename.train  = 'train', filename.pred  = 'pred',
      build.prank     = TRUE, compute.prank = TRUE, write.data.pred = TRUE, write.data.train = TRUE, write.predictor = TRUE)
  }
}  

# random pairing
peptides <- c('ELAGIGILTV', 'GILGFVFTL', 'LLWNGPMAV')

for (peptide in peptides) {
  indir  <- paste0('../step1_', peptide)
  outdir <- file.path(indir, 'TEMPO')
  
  predictor <- TEMPOtrain(input.train = file.path(indir, 'model.csv'), output.path = outdir,
    input.pred = 'validation.csv', filename.train = 'train', filename.pred = 'pred',
    build.prank = TRUE, compute.prank = TRUE, write.data.pred = TRUE, write.data.train = TRUE, write.predictor = TRUE)
}
  
  