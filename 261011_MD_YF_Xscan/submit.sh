#!/bin/bash -l

# submit.sh -- launch TCR-pMHC folding jobs for every method below, all with
# MSAs, over one alanine-scan table.
#
# Set INPUT/METHODS/NUM_SEEDS below, then run: ./submit.sh

METHODS=(AF3 ESMFold2 OpenDDE Protenix)

# One table of TCR-pMHC complexes (cdr3_TRA, cdr3_TRB, TRAV, TRAJ, TRBV, TRBJ,
# peptide, MHC_allele_a, ...), folded as a single job per method -- unlike the
# DuoDesign screens this replaces, which looped over one CSV per target.
# Xscan_ala.csv is 27 complexes: 3 TCRs (TCR_YF, TCR_YF2, TCR_YF3) x 9
# alanine-scan peptides (P1A..P9A, with position 8 already being the WT A).
INPUT=Xscan_ala.csv

# 20 independent seeds, 1 sample per seed. Every method now takes --num-seeds:
# ESMFold2 already did, and the flag was backported to Protenix and OpenDDE
# (their fold-msa.py + *_Launch.sh, from the *_DuoDesign versions that had it;
# originals kept as *.bak-20261011). AF3 takes seeds per *job* rather than per
# run, so for AF3 these go into each inputs/{ID}.json's modelSeeds at build
# time -- AF3_Launch.sh routes the same flag there for you.
NUM_SEEDS=3

# Fixed, rather than each method drawing its own random base seed. Two reasons:
# all four methods then fold the same seed values, and resuming works -- every
# fold-msa.py here decides what to re-run by looking for output named after the
# seeds it was asked for, so a random base seed on a resumed run would match
# nothing and redo the whole batch.
BASE_SEED=1000

# These are the plain TCR pipelines in {method}/, NOT the {method}_DuoDesign/
# ones the previous version of this script used: different input schema (a TCR
# table vs a binder design table + --target) and different MSA source (V/J-gene
# and MHC/B2M MSAs stitched from preMSA/AF3_models, vs one precomputed target
# MSA). AF3/ is new -- it had no TCR counterpart, only AF3_DuoDesign/ -- and
# was written to mirror Protenix/build_input-msa.py exactly, so all four
# methods fold identical sequences from identical MSAs.
METHOD_ROOT=/work/FAC/FBM/LLB/dgfeller/epitope_pred/rroessne

# How each launcher takes its input. Three of them accept
# "<tcr_table.csv> <run_dir>" directly; ESMFold2's historically took a
# pre-built input.csv + msas/ (built by a separate build_input-msa.py run, as
# in the interactive ESMFold2 runs), so a -tcr mode was added to it to build
# those inside the job like the other three. Everything after the run_dir is
# the same for all four.
declare -A LEAD_ARGS=(
    [AF3]=""
    [ESMFold2]="-tcr"
    [OpenDDE]=""
    [Protenix]=""
)

# One wall-clock budget for every method, sized for the slowest. 27 complexes x
# 20 seeds = 540 folds each. Extrapolating the DuoDesign per-design-seed costs
# (ESMFold2 ~5s, AF3 ~17s, Protenix ~17s, OpenDDE ~35s) gives ~45min, ~2.5h,
# ~2.5h and ~5.3h -- but those were smaller complexes than these ~600-residue
# TCR-pMHCs (273 MHC + 108 TCRA + 112 TCRB + 99 B2M + 9 peptide), and each job
# also has to stitch its MSAs up front, so the real numbers will be larger.
# 24h gives OpenDDE several times its estimate; the gpu partition allows
# 3-00:00:00, and jobs are billed for time used rather than time requested, so
# over-asking only risks a slightly longer queue wait. If something does hit
# the limit, just re-run this script: every fold-msa.py here skips complexes a
# previous run already finished and resubmits only the rest.
TIME=10:00:00

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_PATH="$SCRIPT_DIR/$INPUT"

if [ ! -f "$INPUT_PATH" ]; then
    echo "ERROR: input table not found: $INPUT_PATH"
    exit 1
fi

for method in "${METHODS[@]}"; do
    launcher="$METHOD_ROOT/${method}/${method}_Launch.sh"
    if [ ! -f "$launcher" ]; then
        echo "ERROR: no launcher for $method at $launcher"
        exit 1
    fi

    # Each method gets its own run dir, holding input.csv, msas/, inputs/,
    # output/ and output.csv -- and, because we cd into it before sbatch, that
    # job's SLURM log files too.
    run_dir="$SCRIPT_DIR/$method"
    mkdir -p "$run_dir"
    ( cd "$run_dir" || exit 1
      sbatch --time "$TIME" "$launcher" ${LEAD_ARGS[$method]} \
          "$INPUT_PATH" "$run_dir" \
          --num-seeds "$NUM_SEEDS" --base-seed "$BASE_SEED" )
done
