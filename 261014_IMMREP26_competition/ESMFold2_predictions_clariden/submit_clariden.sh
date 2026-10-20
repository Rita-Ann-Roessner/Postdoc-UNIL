#!/bin/bash -l
#
# submit_clariden.sh -- fold the IMMREP26 TCR tables with ESMFold2 on Alps /
# Clariden, one job per epitope. Clariden counterpart of submit.sh (Curnagl).
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
# and benchmark arms were produced. (-msa would route to fold-msa.py.)
#
# Only METHOD_ROOT differs from the Curnagl submit.sh: here it points at the
# Clariden code checkout. The launcher (ESMFold2_Launch.sh) already carries its
# own #SBATCH header for Alps -- partition normal, account prep52, 1 exclusive
# GH200 node (4 GPUs), --environment=esmfold2 -- and splits the CSV over the 4
# GPUs itself, so there is nothing cluster-specific to set here.
#
# Usage:  ./submit_clariden.sh                 # every epitope in EPITOPES below
#         ./submit_clariden.sh GILGFVFTL ...   # just these
#         DRY_RUN=1 ./submit_clariden.sh       # build inputs, submit nothing

set -uo pipefail

# Epitope directories to run. Override by passing names as arguments.
EPITOPES=(CINGVCWTV GILGFVFTL GLCTLVAML LLWNGPMAV LTDEMIAQY NLVPMVATV TTDPSFLGRY YLQPRTFLL)

NUM_SEEDS=3
BASE_SEED=1000

METHOD_ROOT=/iopsstor/scratch/cscs/rroessner/code
BUILD_INPUT="$METHOD_ROOT/ESMFold2/build_input.py"
LAUNCHER="$METHOD_ROOT/ESMFold2/ESMFold2_Launch.sh"

# build_input.py's V/J gene and MHC/B2M reference tables. These live under
# /work/... on Curnagl; on Clariden they are in scratch, so point the script
# there (it reads both from the environment, falling back to the Curnagl paths).
export MIXTCRVIZ_DATA=/iopsstor/scratch/cscs/rroessner/MixTCRviz/data_raw/CDR123
export ALLELES_TAB=/iopsstor/scratch/cscs/rroessner/preMSA/Alleles_tab.txt

# Largest table here is 670 TCRs x 3 seeds = 2010 folds. Billed for time used,
# not requested, and 12 h is the max on the "normal" partition anyway.
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
