#!/bin/bash -l
#
# submit_full.sh -- the FULL-CHECKPOINT arm of the IMMREP26 ESMFold2 runs.
#
# Identical to ../ESMFold2_predictions_clariden (the control) in every respect
# except the checkpoint: biohub/ESMFold2 instead of biohub/ESMFold2-Fast.
# No templates, no pair prompts, no MSAs.
#
# The two checkpoints differ in exactly two config entries; everything else
# (hidden sizes, diffusion module, ESMC-6B backbone, all heads) is the same:
#
#                            Fast      full
#   folding trunk layers       24        48
#   msa_encoder.enabled     False      True
#
# so this arm buys a 2x deeper folding trunk. The MSA encoder is simply skipped
# when no MSA is passed, so enabling it costs nothing here.
#
# Measured: 14.0 s/fold vs 6.8 s on Fast (2.06x). GPU memory is barely changed
# (13.71 GB vs 13.53 GB after load) -- the extra 24 trunk layers are small next
# to ESMC-6B. 7908 folds = ~31 GPU-hours total:
#
#   LTDEMIAQY   390 ->  23 min     NLVPMVATV  1098 ->  64 min
#   GLCTLVAML   576 ->  34 min     TTDPSFLGRY 1161 ->  68 min
#   LLWNGPMAV   687 ->  40 min     YLQPRTFLL  1257 ->  73 min
#   CINGVCWTV   729 ->  42 min     GILGFVFTL  2010 -> 117 min
#
# Seeds are deliberately the SAME as the control arm (base 1000, 3 seeds), so
# the two arms are paired per (ID, seed) and can be compared directly.
#
# Usage:  ./submit_full.sh                 # every epitope in EPITOPES below
#         ./submit_full.sh GILGFVFTL ...   # just these
#         DRY_RUN=1 ./submit_full.sh       # build inputs, submit nothing

set -uo pipefail

EPITOPES=(CINGVCWTV GILGFVFTL GLCTLVAML LLWNGPMAV LTDEMIAQY NLVPMVATV TTDPSFLGRY YLQPRTFLL)

CHECKPOINT=biohub/ESMFold2

NUM_SEEDS=3
BASE_SEED=1000          # same seeds as the control arm -> paired comparison

METHOD_ROOT=/iopsstor/scratch/cscs/rroessner/code
BUILD_INPUT="$METHOD_ROOT/ESMFold2/build_input.py"
LAUNCHER="$METHOD_ROOT/ESMFold2/ESMFold2_Launch.sh"

# build_input.py's V/J gene and MHC/B2M reference tables live under /work on
# Curnagl; on Clariden they are in scratch (it reads both from the environment).
export MIXTCRVIZ_DATA=/iopsstor/scratch/cscs/rroessner/MixTCRviz/data_raw/CDR123
export ALLELES_TAB=/iopsstor/scratch/cscs/rroessner/preMSA/Alleles_tab.txt

# 12 h is the max on "normal"; the largest epitope here is ~2 h.
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
    echo "$epitope: $n complexes x $NUM_SEEDS seeds = $(( n * NUM_SEEDS )) folds" \
         "(~$(( n * NUM_SEEDS * 14 / 4 / 60 )) min on 4 GPUs)"

    if [ -n "$DRY_RUN" ]; then
        echo "   [dry-run] would submit with --checkpoint $CHECKPOINT"
        continue
    fi

    ( cd "$run_dir" && sbatch --time "$TIME" "$LAUNCHER" \
        "$run_dir/input.csv" "$run_dir" \
        --num-seeds "$NUM_SEEDS" --base-seed "$BASE_SEED" \
        --checkpoint "$CHECKPOINT" )
done
