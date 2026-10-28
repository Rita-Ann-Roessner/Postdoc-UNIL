library(TEMPOtrain)

# strings, not numbers: R turns 1.0 into "1" and the filename would not match.
# "all" is the no-selection control: train on every TCR of the epitope.
#FRACS <- c("0.5", "0.8", "1.0", "1.2", "all")
FRACS <- c("all")
EPITOPES <- c("TTDPSFLGRY", "GILGFVFTL", "YLQPRTFLL", "NLVPMVATV", "GLCTLVAML")

TEMPO_DIR <- "TEMPO"                                        # thresholded training sets
PRED_ROOT <- "ESMFold2_predictions_clariden_optimized"      # all TCRs per epitope
OUT_ROOT  <- "TEMPO_out"

for (frac in FRACS) {
  for (ep in EPITOPES) {
    
    pred <- file.path(PRED_ROOT, ep, paste0(ep, ".csv"))
    train <- if (frac == "all") pred else
      file.path(TEMPO_DIR, paste0(ep, "_", frac, ".csv"))
    out <- file.path(OUT_ROOT, paste0(ep, "_", frac))
    
    if (!file.exists(train)) { message("skip, no training set: ", train); next }
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    
    message("=== ", ep, "  frac=", frac, "  (", basename(train), ")")
    TEMPOtrain(input.train = train, output.path = out, input.pred = pred,
               filename.train = ep, filename.pred = "all_tcrs",
               build.prank = TRUE, compute.prank = TRUE,
               write.data.pred = TRUE, write.data.train = TRUE,
               write.predictor = TRUE)
  }
}