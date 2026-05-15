library(TEMPOtrain)
library(MixTCRviz)

peptides <- c("LLWNGPMAV", "ELAGIGILTV", "GILGFVFTL")

for (peptide in peptides) {

  motif_file      <- file.path(peptide, sprintf("A0201_%s.csv", peptide))
  validation_file <- file.path(peptide, "validation.csv")
  output_dir      <- file.path(sprintf("test_TEMPO_validation_%s", peptide))

  TEMPOtrain(
    input.train      = motif_file,
    input.pred       = validation_file,
    output.path      = output_dir,
    filename.train   = "motif_train",
    filename.pred    = "validation_pred",
    build.prank      = FALSE,
    compute.prank    = FALSE,
    write.data.pred  = TRUE,
    write.data.train = FALSE,
    write.predictor  = FALSE
  )

}
