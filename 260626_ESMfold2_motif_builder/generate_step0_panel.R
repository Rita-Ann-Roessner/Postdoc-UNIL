# =============================================================================
# Fixed step-0 panel — generate ONCE, reuse for every epitope.
#
# Produces exactly what run_step0() would produce (same V/J enumeration, same
# flat-baseline CDR3 draw), but stored outside any epitope folder so that every
# epitope folds the IDENTICAL set of TCRs. Because the tcrNNNN ids are then the
# same TCR everywhere, step-0 scores become directly paired across peptides:
# a difference between epitopes is a peptide effect, not a different TCR sample.
#
# ESM_motif_builder.R picks this up via  STEP0_PANEL <- "step0_panel"  and
# re-stamps only the peptide / MHC / species columns per epitope.
#
# Run once:   Rscript generate_step0_panel.R
# Regenerating with a different seed INVALIDATES comparability with everything
# already folded against the old panel — archive, do not overwrite.
# =============================================================================

.sourced_for_benchmark <- TRUE
source("ESM_motif_builder.R")   # functions + setup (cdr3_baseline, INPUT_DIR); main skipped

OUT_DIR <- "step0_panel"
SEED    <- 20260101             # fixed and recorded; the panel's identity

if (dir.exists(OUT_DIR) && length(list.files(OUT_DIR, pattern = "^model_.*\\.csv$"))) {
  stop(sprintf(paste0("'%s' already contains a panel.\n",
                      "Overwriting would silently break comparability with every\n",
                      "epitope already folded against it. Archive it first if you\n",
                      "really mean to replace it."), OUT_DIR))
}
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(SEED)

# Placeholder peptide/MHC: these columns are overwritten per epitope by run_step0,
# so their values here are irrelevant — only the TCR columns and ids are reused.
run_step0(
  peptide         = "XXXXXXXXX",
  mhc_allele      = "HLA_A0201",
  label           = ".",                 # writes straight into OUT_DIR/./step0
  cdr3_baseline   = cdr3_baseline,
  base_output_dir = OUT_DIR,
  input_dir       = INPUT_DIR,
  panel_dir       = NULL                 # generate, do not read a panel
)

# lift the two model tables out of the nested step0/ dir into OUT_DIR, then drop
# the scratch dir: its *_seqs.csv carry the placeholder peptide and must never be
# folded by accident.
for (cn in c("alpha", "beta")) {
  from <- file.path(OUT_DIR, "step0", sprintf("model_%s.csv", cn))
  file.copy(from, file.path(OUT_DIR, sprintf("model_%s.csv", cn)), overwrite = TRUE)
}
unlink(file.path(OUT_DIR, "step0"), recursive = TRUE)

writeLines(c(sprintf("seed: %d", SEED),
             sprintf("created: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
             sprintf("n_alpha: %d", nrow(read.csv(file.path(OUT_DIR, "model_alpha.csv")))),
             sprintf("n_beta: %d",  nrow(read.csv(file.path(OUT_DIR, "model_beta.csv")))),
             sprintf("CDR3_TEMPLATED at creation: %s", CDR3_TEMPLATED),
             sprintf("JUNCTION_PSSM at creation: %s", JUNCTION_PSSM),
             "cdr3 generation: draw_random_cdr3() -> assemble_templated_cdr3()",
             "  (germline V block + junction + germline J block), the same path",
             "  the enrichment steps use."),
           file.path(OUT_DIR, "PANEL_INFO.txt"))

message(sprintf("\nPanel written to %s/  (model_alpha.csv, model_beta.csv, PANEL_INFO.txt)", OUT_DIR))
message("ESM_motif_builder.R will use it for every epitope via STEP0_PANEL.")
