#!/bin/bash -l

# submit.sh -- launch DuoDesign folding jobs for every method x target
# combination below.
#
# Set TARGETS/METHODS below, then run: ./submit.sh

METHODS=(AF3 ESMFold2 OpenDDE Protenix)

TARGETS=(IL-6_WT IL-6_F102A IL-6_E138A IL-6_F153A)

# 20 independent seeds, 1 sample per seed -- forwarded to every method's
# fold-msa.py (and, for AF3, also to build_input-msa.py's modelSeeds) via the
# shared --num-seeds flag now supported by all four *_Launch.sh scripts.
NUM_SEEDS=20

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

# One wall-clock budget for every method, sized for the slowest. Only the
# folding scales with NUM_SEEDS, so this is extrapolated from the 3-seed run
# (jobs 65513757-72) by fitting elapsed time against each target's design
# count, which separates the fixed MSA/model-load overhead from the per-fold
# cost:
#
#   per design-seed: ESMFold2 ~5s, AF3 ~17s, Protenix ~17s, OpenDDE ~35s
#
# The largest job is IL-6_WT at 17 designs, i.e. 17*20 = 340 folds, giving
# ~0.5h (ESMFold2), ~1.6h (AF3, Protenix) and ~3.3h (OpenDDE). 8h covers
# OpenDDE with ~2.4x margin and is simply generous for the rest; the gpu
# partition allows 3-00:00:00, and jobs are billed for time used rather than
# time requested, so over-asking only risks a slightly longer queue wait.
TIME=08:00:00

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for method in "${METHODS[@]}"; do
    ( cd "$SCRIPT_DIR/${LOCAL_DIR[$method]}" || exit 1
      for target in "${TARGETS[@]}"; do
          sbatch --time "$TIME" "$DUODESIGN_ROOT/${method}_DuoDesign/${method}_Launch.sh" \
              "../${target}.csv" "./${target}" --target "${target}" --num-seeds "${NUM_SEEDS}"
      done )
done
