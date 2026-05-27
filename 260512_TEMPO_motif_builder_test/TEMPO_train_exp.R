library(TEMPOtrain)
library(MixTCRviz)
library(pROC)

# =============================================================================
# TEMPO upper-ceiling benchmark
# Trains TEMPO directly on known binders and scores the validation set.
# AUC / AUC0.1 are computed per epitope and saved to auc_summary_exp.csv.
# =============================================================================

EPITOPES_FILE <- "epitopes.txt"
epitopes <- trimws(readLines(EPITOPES_FILE))
epitopes <- epitopes[nchar(epitopes) > 0]

auc_summary <- data.frame()

for (epitope in epitopes) {

  known_binders_file <- file.path(epitope, paste0(epitope, ".csv"))
  validation_file    <- file.path(epitope, "validation.csv")
  output_dir         <- epitope

  if (!file.exists(known_binders_file)) {
    message(sprintf("Skipping %s — known binders file not found", epitope))
    next
  }
  if (!file.exists(validation_file)) {
    message(sprintf("Skipping %s — validation file not found", epitope))
    next
  }

  message(sprintf("Processing %s ...", epitope))

  tryCatch({
    TEMPOtrain(
      input.train      = known_binders_file,
      input.pred       = validation_file,
      output.path      = output_dir,
      filename.pred    = "TEMPO_train_exp",
      build.prank      = FALSE,
      compute.prank    = FALSE,
      write.data.pred  = TRUE,
      write.data.train = FALSE,
      write.predictor  = FALSE
    )

    pred_files <- list.files(output_dir, pattern = "^TEMPO_train_exp", full.names = TRUE)
    if (length(pred_files) == 0) {
      message(sprintf("  No prediction output found for %s — skipping AUC", epitope))
      next
    }
    pred_df <- read.csv(pred_files[1])

    if (!"Label" %in% colnames(pred_df) || !"score" %in% colnames(pred_df)) {
      message(sprintf("  Missing Label or score column for %s — skipping AUC", epitope))
      next
    }

    pred_df <- pred_df[!is.na(pred_df$score), ]
    roc_obj   <- roc(pred_df$Label, pred_df$score, quiet = TRUE)
    auc_val   <- as.numeric(auc(roc_obj))
    auc01_val <- as.numeric(auc(roc_obj,
                                partial.auc         = c(1, 0.9),
                                partial.auc.correct = TRUE,
                                partial.auc.focus   = "specificity"))

    message(sprintf("  AUC = %.4f  |  AUC0.1 = %.4f", auc_val, auc01_val))

    auc_summary <- rbind(auc_summary, data.frame(
      epitope = epitope,
      auc     = auc_val,
      auc01   = auc01_val,
      stringsAsFactors = FALSE
    ))

  }, error = function(e) {
    message(sprintf("  ERROR for %s: %s", epitope, e$message))
  })
}

write.csv(auc_summary, "auc_summary_exp.csv", row.names = FALSE)

message("\n===== AUC Summary (experimental) =====")
for (i in seq_len(nrow(auc_summary))) {
  message(sprintf("  %s: AUC = %.4f  |  AUC0.1 = %.4f",
                  auc_summary$epitope[i], auc_summary$auc[i], auc_summary$auc01[i]))
}
message(sprintf("  Saved to: auc_summary_exp.csv"))
