#!/bin/bash -l
#
# submit.sh -- fold the IMMREP26 TCR tables with ESMFold2, one job per epitope.
#
# Each <epitope>/<epitope>.csv is a raw TCR table (cdr3_TRA, cdr3_TRB, TRAV,
# TRAJ, TRBV, TRBJ, peptide, MHC_allele_a). fold.py needs sequences instead
# (ID, TCRA, TCRB, MHC, B2M, PEPTIDE), so build_input.py writes <epitope>/
# input.csv here, at submit time -- it is plain stdlib python3 and takes
# seconds, so a bad gene name fails now rather than after a queue wait on a
# GPU.
#
# Everything for an epitope lives in its own directory: the raw table, the
# built input.csv, models/, metrics/, output.csv and the job's SLURM logs.
#
# No MSAs: the bare launcher form runs fold.py, matching how the ESMFold2 atlas
# and benchmark arms were produced. (-tcr/-msa would route to fold-msa.py.)
#
# Usage:  ./submit.sh                 # every epitope in EPITOPES below
#         ./submit.sh GILGFVFTL ...   # just these
#         DRY_RUN=1 ./submit.sh       # build inputs, submit nothing

set -uo pipefail

# Epitope directories to run. Override by passing names as arguments.
#EPITOPES=(GILGFVFTL GLCTLVAML NLVPMVATV TTDPSFLGRY YLQPRTFLL)
EPITOPES=(LTDEMIAQY CINGVCWTV LLWNGPMAV)

NUM_SEEDS=3
BASE_SEED=1000

METHOD_ROOT=/work/FAC/FBM/LLB/dgfeller/epitope_pred/rroessne
BUILD_INPUT="$METHOD_ROOT/ESMFold2/build_input.py"
LAUNCHER="$METHOD_ROOT/ESMFold2/ESMFold2_Launch.sh"

# Largest table here is 670 TCRs x 3 seeds = 2010 folds. Billed for time used,
# not requested, so this is deliberately generous.
TIME=12:00:00

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DRY_RUN="${DRY_RUN:-}"

[ $# -gt 0 ] && EPITOPES=("$@")

for epitope in "${EPITOPES[@]}"; do
    run_dir="$SCRIPT_DIR/$epitope"
    src="$run_dir/$epitope.csv"
    if [ ! -f "$src" ]; then
        echo "[skip] $epitope: no $epitope.csv"
        continue
    fi

    if ! python3 "$BUILD_INPUT" "$src" "$run_dir/input.csv"; then
        echo "[skip] $epitope: build_input.py failed"
        continue
    fi

    n=$(( $(wc -l < "$run_dir/input.csv") - 1 ))
    echo "$epitope: $n complexes x $NUM_SEEDS seeds = $(( n * NUM_SEEDS )) folds"

    if [ -n "$DRY_RUN" ]; then
        echo "   [dry-run] would submit"
        continue
    fi

    ( cd "$run_dir" && sbatch --time "$TIME" "$LAUNCHER" \
        "$run_dir/input.csv" "$run_dir" \
        --num-seeds "$NUM_SEEDS" --base-seed "$BASE_SEED" )
done
