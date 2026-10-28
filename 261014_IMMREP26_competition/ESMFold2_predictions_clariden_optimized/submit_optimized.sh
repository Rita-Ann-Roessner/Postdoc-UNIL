#!/bin/bash -l
#
# submit_optimized.sh -- ESMFold2_predictions_clariden_optimized
#
# Same complexes and seeds as every other arm; the differences are:
#
#   num_loops           3  -> 20    (the model's own config default)
#   num_sampling_steps 50  -> 200   (ditto)
#   checkpoint                      biohub/ESMFold2-Fast
#   + ipSAE scores and the raw PAE matrix are written
#
# WHY THESE TWO KNOBS
#   num_loops is recycling: the trunk runs num_loops+1 times, each pass
#   re-reading its own pair state. It improves WHAT the model believes, which
#   is the lever for docking -- the part we found to be the bottleneck when
#   pinning the pMHC to 0.18 Ang did not move discrimination.
#   num_sampling_steps is diffusion denoising: it improves how precisely that
#   belief is rendered into coordinates, not the belief itself.
#
# MEASURED COST (biohub/ESMFold2-Fast at 20/200): 35 s/fold, vs 6.8 s at 3/50.
#
#   LTDEMIAQY     390 folds ->  0.9 h on 4 GPUs
#   GLCTLVAML     576 folds ->  1.4 h on 4 GPUs
#   LLWNGPMAV     687 folds ->  1.7 h on 4 GPUs
#   CINGVCWTV     729 folds ->  1.8 h on 4 GPUs
#   NLVPMVATV    1098 folds ->  2.7 h on 4 GPUs
#   TTDPSFLGRY   1161 folds ->  2.8 h on 4 GPUs
#   YLQPRTFLL    1257 folds ->  3.1 h on 4 GPUs
#   GILGFVFTL    2010 folds ->  4.9 h on 4 GPUs
#   total 7908 folds ~ 77 GPU-hours
#
# The 12 h limit on "normal" is the binding constraint for GILGFVFTL; if a node
# is slow, split it with code/ESMFold2/split_and_submit.sh.
#
# ipSAE (Dunbrack 2025) is written alongside ipTM: an asymmetric per-chain-pair
# matrix plus ipsae_pair_mean, using the d0dom variant, and the merged
# TCR(A+B) -> MHC+PEPTIDE pairings. Raw PAE goes to pae/<ID>.npz (~0.6 MB each,
# ~5 GB per arm) so other cutoffs or scores can be computed without refolding.
#
# Usage:  ./submit_optimized.sh                 # every epitope
#         ./submit_optimized.sh GILGFVFTL ...   # just these
#         DRY_RUN=1 ./submit_optimized.sh

set -uo pipefail

EPITOPES=(CINGVCWTV GILGFVFTL GLCTLVAML LLWNGPMAV LTDEMIAQY NLVPMVATV TTDPSFLGRY YLQPRTFLL)

CHECKPOINT=biohub/ESMFold2-Fast
NUM_LOOPS=20
NUM_SAMPLING_STEPS=200

NUM_SEEDS=3
BASE_SEED=1000          # same seeds as every other arm -> paired comparison

METHOD_ROOT=/iopsstor/scratch/cscs/rroessner/code
BUILD_INPUT="$METHOD_ROOT/ESMFold2/build_input.py"
LAUNCHER="$METHOD_ROOT/ESMFold2/ESMFold2_Launch.sh"

export MIXTCRVIZ_DATA=/iopsstor/scratch/cscs/rroessner/MixTCRviz/data_raw/CDR123
export ALLELES_TAB=/iopsstor/scratch/cscs/rroessner/preMSA/Alleles_tab.txt

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
         "(~$(( n * NUM_SEEDS * 35 / 4 / 3600 )) h on 4 GPUs)"

    if [ -n "$DRY_RUN" ]; then
        echo "   [dry-run] would submit --checkpoint $CHECKPOINT --num-loops $NUM_LOOPS" \
             "--num-sampling-steps $NUM_SAMPLING_STEPS --save-pae"
        continue
    fi

    ( cd "$run_dir" && sbatch --time "$TIME" "$LAUNCHER" \
        "$run_dir/input.csv" "$run_dir" \
        --num-seeds "$NUM_SEEDS" --base-seed "$BASE_SEED" \
        --checkpoint "$CHECKPOINT" \
        --num-loops "$NUM_LOOPS" --num-sampling-steps "$NUM_SAMPLING_STEPS" \
        --save-pae )
done
