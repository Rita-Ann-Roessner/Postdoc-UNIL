devtools::load_all("../TEMPOtrain") # treat the TEMPOtrain folder like an R package and load it
library(TEMPOtrain)

# train on celiac specific TCRs
epitopes <- c('FPQPEQPFPWQP', 'PQPELPYPQPE', 'PQQPFPQPEQPFP', 'QLQPFPQPELPY')

for (epitope in epitopes){
  
  predictor <- TEMPOtrain(input.train = paste0("celiac_specific_CD4_", epitope,".csv"), output.path = paste0("TEMPO_out/", epitope),
              input.pred = "celiac_CD4_CD8_clean.csv", filename.train = "celiac_specific", filename.pred="repertoire", build.prank = T, compute.prank=T,
              write.data.pred = T, write.data.train = T, write.predictor = T)
}


# apparently I can save the predictor (with the built prank and then lod it into input.pred=predictor$predictor)