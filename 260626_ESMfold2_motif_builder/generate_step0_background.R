# =============================================================================
# Step-0 background (decoy) repertoire
#
# Folds the FIXED step-0 panel (step0_panel/, via STEP0_PANEL) against a decoy
# peptide. Because it is the same panel every cognate epitope uses, each TCR is
# scored against both its cognate peptide and the decoy, so
#       score(t, cognate) - score(t, decoy)
# is a paired, per-TCR contrast with no CDR3 sampling noise in it at all.
#
# No replicates: replicates existed to average over independent CDR3 draws, which
# the fixed panel makes unnecessary for a paired contrast (and the panel already
# carries ~50 alpha TCRs per V gene). If you ever do need CDR3-averaged per-gene
# rates -- mainly for beta, which only gets ~13 TCRs per gene -- set
# STEP0_PANEL <- NULL and loop this script over several seeds instead.
#
# Decoy = a neutral, non-cognate peptide (A2 anchors + poly-Ala TCR face,
# "ALAAAAAAV"), so the background isolates the epitope-INDEPENDENT fold artifact
# rather than real germline-encoded specificity. Works for any peptide, though.
#
# NOTE: these folds inherit the current CDR3_TEMPLATED setting. Do not pool them
# with the archived pre-templating replicates in bck.step0_background/.
#
# Workflow:
#   Rscript generate_step0_background.R        # write the decoy batch
#   (fold <OUT_DIR>/<label>/step0/model_{alpha,beta}_seqs.csv on the cluster;
#    place results as output_{alpha,beta}.csv beside the models)
#   Rscript analyze_step0_background.R         # -> bg(g) tables + artifact check
# =============================================================================

.sourced_for_benchmark <- TRUE
source("ESM_motif_builder.R")   # functions + setup (cdr3_baseline, INPUT_DIR); main skipped

# ---- config -----------------------------------------------------------------
OUT_DIR     <- "step0_background"
# "MHC_PEPTIDE" labels (same convention as `epitopes`). Decoy peptide:
#   ALAAAAAAV  neutral, featureless poly-Ala TCR face (A2 anchors L2/V9)
BG_PEPTIDES <- c("A0201_ALAAAAAAV")

if (is.null(STEP0_PANEL) || !nzchar(STEP0_PANEL))
  stop("STEP0_PANEL is not set. This script folds the fixed panel against the decoy;\n",
       "without a panel the background would not be paired with the cognate runs.")
if (!file.exists(file.path(STEP0_PANEL, "model_alpha.csv")))
  stop(sprintf("Panel '%s' not found. Create it once with: Rscript generate_step0_panel.R",
               STEP0_PANEL))

# ---- generate ---------------------------------------------------------------
for (ep in BG_PEPTIDES) {
  mhc        <- sub("_.*", "", ep)
  peptide    <- sub("^[^_]*_", "", ep)
  mhc_allele <- if (grepl("^HLA_", mhc)) mhc else paste0("HLA_", mhc)

  message(sprintf("\n[%s] decoy background from panel '%s'", ep, STEP0_PANEL))
  run_step0(
    peptide         = peptide,
    mhc_allele      = mhc_allele,
    label           = ep,
    cdr3_baseline   = cdr3_baseline,
    base_output_dir = OUT_DIR,
    input_dir       = INPUT_DIR,
    panel_dir       = STEP0_PANEL      # explicit: same TCRs as every cognate run
  )
}

message(sprintf(
  paste0("\nDone. Fold on the cluster:\n  %s/<label>/step0/model_{alpha,beta}_seqs.csv\n",
         "then place results as output_{alpha,beta}.csv beside the models and run:\n",
         "  Rscript analyze_step0_background.R\n\n",
         "DECOY_DIR in ESM_motif_builder.R must point at '%s'."),
  OUT_DIR, OUT_DIR))
