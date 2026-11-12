#!/bin/bash -l

# submit_bbf14.sh -- launch BBF-14 DuoDesign folding jobs, all 4 methods,
# with no target MSA (folds BBF-14 single-sequence via each method's
# --no-msa). Independent of submit.sh's per-METHOD/per-target loop; run this
# whenever you want the no-MSA BBF-14 comparison.

TARGET=BBF-14

# 3 independent seeds, 1 sample per seed -- forwarded to every method's
# fold-msa.py (and, for AF3, also to build_input-msa.py's modelSeeds) via the
# shared --num-seeds flag supported by all four *_Launch.sh scripts.
NUM_SEEDS=3

DUODESIGN_ROOT=/work/FAC/FBM/LLB/dgfeller/epitope_pred/rroessne

# Local subdirectory (relative to this script) to cd into before sbatch, so
# ../BBF-14.csv and ./BBF-14_noMSA resolve correctly. Matches submit.sh's
# LOCAL_DIR mapping (Protenix_v2 is Protenix_Launch.sh's default model since
# no model_name is passed below).
declare -A LOCAL_DIR=(
    [AF3]=AF3
    [ESMFold2]=ESMFold2
    [OpenDDE]=OpenDDE
    [Protenix]=Protenix_v2
)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#for m in AF3 ESMFold2 OpenDDE Protenix; do
for m in AF3; do
    (
        cd "$SCRIPT_DIR/${LOCAL_DIR[$m]}" || exit 1
        # --time bumped 3x over the single-seed 00:30:00 budget since each
        # job now folds NUM_SEEDS times as much -- re-check against actual
        # single-seed wall time before relying on this for a bigger batch.
        sbatch --time 01:30:00 "$DUODESIGN_ROOT/${m}_DuoDesign/${m}_Launch.sh" \
            "../${TARGET}.csv" "./${TARGET}_noMSA" --target "${TARGET}" --no-msa --num-seeds "${NUM_SEEDS}"
    )
done
