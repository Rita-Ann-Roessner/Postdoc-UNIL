# =============================================================================
# TRAV17 / TRAJ31 test panel
#
# Question: in the NY-ESO-1 step-0 sweep, the single TRAV17/TRAJ31 fold scored
# 0.818 -- rank 8 of 2267 V/J cells -- while TRAV17 averaged over all its TRAJ
# ranked only 26 of 45 genes. One fold cannot distinguish
#
#   (a) TRAV17/TRAJ31 is a genuinely strong pair that the factorized
#       P(v) x P(j) draw almost never visits, from
#   (b) it was the top of 51 single draws, i.e. noise.
#
# Design: two arms, TRAV17 held fixed, dummy beta partner as in the pipeline.
#   arm "TRAJ31" : J fixed to TRAJ31
#   arm "otherJ" : J drawn uniformly from the other TRAV17 cells that have a
#                  CDR3 baseline (mirrors the step-0 sweep, where every cell
#                  appears exactly once, so J is effectively uniform)
# Both arms are matched on CDR3 length, because ipTM falls with alpha CDR3
# length (Spearman -0.28 over the step-0 panel) and the real binders are short.
# Without the match, a TRAJ31 effect and a length effect are not separable.
#
# CDR3s come from draw_random_cdr3() in ESM_motif_builder.R, UNMODIFIED, so the
# germline-templated construction is identical to step 0 and the enrichment
# steps. Length is controlled by rejection: draw, keep if the length is one we
# want. Nothing about the sampling path is changed, only which draws are kept.
#
# Usage:  Rscript generate_TRAJ31_test.R
# =============================================================================

.sourced_for_benchmark <- TRUE          # skip the cluster-in-the-loop main run
source("ESM_motif_builder.R")           # brings in draw_random_cdr3, cdr3_baseline,
                                        # GENE_TABLES, ALLELES_SEQ, DUMMY_TRB,
                                        # write_fold_input_csv, CDR3_TEMPLATED

### ---- Configuration --------------------------------------------------------

RUN_DIR      <- "TCR_motif_atlas_cdr3_templated_new_installation_reduced_cycles_NY-ESO-1"
OUT_DIR      <- file.path(RUN_DIR, "struc_preds_traj31_test")
V_GENE       <- "TRAV17"
TEST_J       <- "TRAJ31"
LENGTHS      <- c(11, 12, 13, 14)   # 11 = the real binders; 13 = the 0.818 fold
N_PER_CELL   <- 10                  # folds per arm per length -> 40 per arm, 80 total
PEPTIDE      <- "SLLMWITQC"
MHC_ALLELE   <- "HLA_A0201"
SPECIES      <- "HomoSapiens"
DEDUP        <- TRUE                # drop repeat CDR3s within an arm x length cell:
                                    # folding the same sequence twice yields nothing.
                                    # Set FALSE to sample strictly as the pipeline does.
MAX_TRIES    <- 20000               # rejection-sampling budget per cell
SEED         <- 1

set.seed(SEED)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

### ---- Which J genes are eligible -------------------------------------------
# Same gate the generator uses: the (V, J) pair must have a CDR3 baseline.

all_keys  <- names(cdr3_baseline[["TRA"]])
v_keys    <- grep(paste0("^", V_GENE, "_"), all_keys, value = TRUE)
js_all    <- sub(paste0("^", V_GENE, "_"), "", v_keys)
js_other  <- setdiff(js_all, TEST_J)

if (!(TEST_J %in% js_all))
  stop(sprintf("%s_%s has no CDR3 baseline.", V_GENE, TEST_J))

message(sprintf("%s: %d V/J cells with a CDR3 baseline (%d available to the control arm)",
                V_GENE, length(js_all), length(js_other)))

### ---- Draw CDR3s at a target length ----------------------------------------
# draw_random_cdr3() samples the length itself from the V/J/L baseline, so we
# draw repeatedly and keep the hits. The sampling path is untouched.

draw_at_length <- function(j_pool, target_len, n_want) {
  out <- list(); tries <- 0L
  while (length(out) < n_want && tries < MAX_TRIES) {
    tries <- tries + 1L
    j   <- if (length(j_pool) == 1L) j_pool else sample(j_pool, 1L)
    s   <- draw_random_cdr3("TRA", V_GENE, j, cdr3_baseline)
    if (is.na(s) || nchar(s) != target_len) next
    if (DEDUP && any(vapply(out, function(o) identical(o$cdr3, s), logical(1)))) next
    out[[length(out) + 1L]] <- list(j = j, cdr3 = s)
  }
  if (length(out) < n_want)
    warning(sprintf("  L=%d: only %d/%d drawn after %d tries", target_len, length(out), n_want, tries))
  out
}

### ---- Build both arms ------------------------------------------------------

rows <- list()
for (arm in c("TRAJ31", "otherJ")) {
  pool <- if (arm == "TRAJ31") TEST_J else js_other
  for (L in LENGTHS) {
    got <- draw_at_length(pool, L, N_PER_CELL)
    message(sprintf("  arm %-7s L=%2d -> %d CDR3s", arm, L, length(got)))
    for (g in got) {
      rows[[length(rows) + 1L]] <- data.frame(
        arm          = arm,
        TRAV         = V_GENE,
        TRAJ         = g$j,
        cdr3_TRA     = g$cdr3,
        cdr3_len     = nchar(g$cdr3),
        TRBV         = DUMMY_TRB$v,
        TRBJ         = DUMMY_TRB$j,
        cdr3_TRB     = DUMMY_TRB$cdr3,
        peptide      = PEPTIDE,
        MHC_allele_a = MHC_ALLELE,
        species      = SPECIES,
        stringsAsFactors = FALSE)
    }
  }
}
panel <- do.call(rbind, rows)
panel$id <- sprintf("tcr_%04d", seq_len(nrow(panel)))

### ---- Write fold input + annotation ----------------------------------------

name     <- "traj31_test"
fold_csv <- file.path(OUT_DIR, paste0(name, ".csv"))
write_fold_input_csv(panel, fold_csv)

ann <- panel[, c("id", "arm", "TRAV", "TRAJ", "cdr3_TRA", "cdr3_len",
                 "TRBV", "TRBJ", "cdr3_TRB")]
names(ann)[1] <- "ID"
write.csv(ann, file.path(OUT_DIR, paste0(name, "_annotation.csv")), row.names = FALSE)

### ---- submit.sh -------------------------------------------------------------
# Same shape as struc_preds_dummy_beta/submit.sh: the CSV already holds complete
# sequences, so it goes straight to fold.py with the reduced-cycle settings the
# rest of this run used (num_loops=3, num_sampling_steps=50). Point SCRATCH at
# wherever you stage it on the cluster.

submit <- c(
  "#!/bin/bash",
  "# ESMFold2-Fast predictions for the TRAV17/TRAJ31 test panel, reduced cycles",
  "# (num_loops=3, num_sampling_steps=50) -- MUST match the step-0 settings, or",
  "# the scores are not comparable to the panel this test is interpreted against.",
  "",
  "RUN_DIR=/iopsstor/scratch/cscs/rroessner/struc_preds_traj31_test",
  "LAUNCHER=/iopsstor/scratch/cscs/rroessner/code/ESMFold2/ESMFold2_Launch.sh",
  "",
  "NUM_LOOPS=3",
  "NUM_SAMPLING_STEPS=50",
  "",
  sprintf("CSVS=(\n    %s.csv\n)", name),
  "",
  "cd \"$RUN_DIR\" || exit 1",
  "",
  "for csv in \"${CSVS[@]}\"; do",
  "    if [ ! -f \"$csv\" ]; then",
  "        echo \"[skip] $csv not found\"",
  "        continue",
  "    fi",
  "    n=\"${csv%.csv}\"",
  "    echo \"[submit] $n\"",
  "    sbatch \"$LAUNCHER\" \"$RUN_DIR/$csv\" \"$RUN_DIR/$n\" \\",
  "        --num-loops \"$NUM_LOOPS\" --num-sampling-steps \"$NUM_SAMPLING_STEPS\"",
  "done")
writeLines(submit, file.path(OUT_DIR, "submit.sh"))
Sys.chmod(file.path(OUT_DIR, "submit.sh"), "0755")

### ---- Summary --------------------------------------------------------------

message("\n=== panel ===")
print(table(panel$arm, panel$cdr3_len))
message(sprintf("\ndistinct CDR3s: %d of %d rows", length(unique(panel$cdr3_TRA)), nrow(panel)))
message(sprintf("distinct J in control arm: %d", length(unique(panel$TRAJ[panel$arm == "otherJ"]))))
message(sprintf("\nwrote %s\n      %s", fold_csv,
                file.path(OUT_DIR, paste0(name, "_annotation.csv"))))
