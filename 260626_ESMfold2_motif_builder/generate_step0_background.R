# =============================================================================
# Step-0 background replicates
#
# Generate N independent flat-baseline step-0 repertoires (identical V/J grid,
# FRESH random CDR3s each replicate) for a decoy peptide. Fold all replicates on
# the cluster, then average per-V / per-length ESMFold pass rates across them
# (analyze_step0_background.R) to get a low-noise, CDR3-agnostic background bg(g).
#
# Decoy = a neutral, non-cognate peptide (Option C: A2 anchors + poly-Ala TCR
# face, "ALAAAAAAV"), so bg(g) isolates the epitope-INDEPENDENT ESMFold fold
# artifact (what we want to subtract) rather than real germline-encoded
# specificity for a cognate epitope. Works for any peptide, though.
#
# Workflow:
#   Rscript generate_step0_background.R        # make rep1..N step-0 batches
#   (fold each rep's model_{alpha,beta}_seqs.csv on the cluster; place results
#    as output_{alpha,beta}.csv beside the models)
#   Rscript analyze_step0_background.R         # pool -> bg(g) table + artifact check
# =============================================================================

.sourced_for_benchmark <- TRUE
source("ESM_motif_builder.R")   # functions + setup (cdr3_baseline, INPUT_DIR); main skipped

# ---- config -----------------------------------------------------------------
OUT_DIR     <- "step0_background"
N_REPS      <- 3
BASE_SEED   <- 100
# "MHC_PEPTIDE" labels (same convention as `epitopes`). Decoy peptide:
#   ALAAAAAAV  Option C: neutral, featureless poly-Ala TCR face (anchors L2/V9)
BG_PEPTIDES <- c("A0201_ALAAAAAAV")

# ---- generate ---------------------------------------------------------------
for (r in seq_len(N_REPS)) {
  set.seed(BASE_SEED + r)                       # distinct, reproducible repertoire
  rep_dir <- file.path(OUT_DIR, sprintf("rep%d", r))

  for (ep in BG_PEPTIDES) {
    mhc        <- sub("_.*", "", ep)
    peptide    <- sub("^[^_]*_", "", ep)
    mhc_allele <- if (grepl("^HLA_", mhc)) mhc else paste0("HLA_", mhc)

    message(sprintf("\n[rep %d/%d | %s] generating step-0 background repertoire",
                    r, N_REPS, ep))
    run_step0(
      peptide         = peptide,
      mhc_allele      = mhc_allele,
      label           = ep,
      cdr3_baseline   = cdr3_baseline,
      base_output_dir = rep_dir,
      input_dir       = INPUT_DIR
    )
  }
}

message(sprintf(
  "\nDone. Fold each replicate on the cluster:\n  %s/rep*/<label>/step0/model_{alpha,beta}_seqs.csv\nthen place ESMFold results as output_{alpha,beta}.csv beside the models and run:\n  Rscript analyze_step0_background.R",
  OUT_DIR))
