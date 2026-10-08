#!/bin/bash -l

# submit.sh -- launch DuoDesign folding jobs for every method x target
# combination below.
#
# Set TARGETS/METHODS below, then run: ./submit.sh

METHODS=(AF3 ESMFold2 OpenDDE Protenix)

TARGETS=(IL-6_WT IL-6_F102A IL-6_E138A IL-6_F153A)

# 3 independent seeds, 1 sample per seed -- forwarded to every method's
# fold-msa.py (and, for AF3, also to build_input-msa.py's modelSeeds) via the
# shared --num-seeds flag now supported by all four *_Launch.sh scripts.
NUM_SEEDS=3

DUODESIGN_ROOT=/work/FAC/FBM/LLB/dgfeller/epitope_pred/rroessne

# Local subdirectory (relative to this script) to cd into before sbatch, so
# ../<target>.csv and ./h<target> resolve correctly. Matches Protenix_Launch.sh's
# default model (protenix-v2) since no model_name is passed below.
declare -A LOCAL_DIR=(
    [AF3]=AF3
    [ESMFold2]=ESMFold2
    [OpenDDE]=OpenDDE
    [Protenix]=Protenix_v2
)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for method in "${METHODS[@]}"; do
    ( cd "$SCRIPT_DIR/${LOCAL_DIR[$method]}" || exit 1
      for target in "${TARGETS[@]}"; do
          # --time bumped 3x over the original 00:30:00 single-seed budget
          # since each job now folds NUM_SEEDS times as much -- re-check
          # against actual single-seed wall time before relying on this for
          # a bigger batch.
          sbatch --time 01:30:00 "$DUODESIGN_ROOT/${method}_DuoDesign/${method}_Launch.sh" \
              "../${target}.csv" "./${target}" --target "${target}" --num-seeds "${NUM_SEEDS}"
      done )
done
