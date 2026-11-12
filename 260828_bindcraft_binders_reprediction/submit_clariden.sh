#!/bin/bash -l

# submit.sh -- ESMFold2 repredictions of the BindCraft binders/decoys, on
# Alps/Clariden. One script, one sbatch job per target, all four targets:
#
#   DerF7, DerF21, IFNAR2   with the target's precomputed MSA
#   BBF-14                  single-sequence (--no-msa)
#
# BBF-14 is the odd one out on purpose: it is the no-evolutionary-information
# control, so its target chain is folded single-sequence. The MSA is the ONLY
# difference -- same checkpoint, same num_loops/num_sampling_steps, same
# seeds -- so the four runs stay directly comparable. (The Curnagl original
# split this across submit.sh + submit_bbf14_template.sh because the other
# three predictors needed a manually supplied 9HAC template for the no-MSA
# arm; ESMFold2 has no template channel at all, so there is nothing to
# configure and the two scripts collapse into this one.)
#
# Only ESMFold2 is being re-run here -- AF3/OpenDDE/Protenix are not ported to
# Clariden in this pipeline, so there is no METHODS loop any more.
#
# Each job gets one exclusive Clariden node (4x GH200) and ESMFold2_Launch.sh
# splits that target's designs over the 4 GPUs. Run: ./submit.sh

# Targets folded WITH the precomputed target MSA (preMSA/<Target>.a3m in the
# code dir). Output dir per target: ESMFold2/<Target>/
MSA_TARGETS=(IFNAR2 DerF7 DerF21)

# Targets folded WITHOUT any target MSA. Output dir: ESMFold2/<Target>_noMSA/
NOMSA_TARGETS=(BBF-14)

# 3 independent seeds, 1 sample per seed -- forwarded to fold-msa.py.
NUM_SEEDS=3

# fold-msa.py's defaults are the model config's own values (num_loops 20,
# num_sampling_steps 200), so nothing is passed here. The Curnagl copy of
# fold-msa.py hardcoded 3/50 instead; anything predicted with that is not
# comparable with these runs.
EXTRA_ARGS=()

# Wall time per job. Measured on Clariden 2026-11-04, biohub/ESMFold2 at
# num_loops=20/num_sampling_steps=200: 15 s per fold for the largest complex
# here (DerF21 target 118 aa + 183 aa binder = 301 residues; 15.3 GB peak GPU,
# so 4 folds share a GH200 node comfortably), plus ~30 s of model load per GPU.
# Worst target is BBF-14 at 10 designs: split over 4 GPUs that is 3 designs x
# 3 seeds x 15 s + load, under 3 min. 30 min is therefore a ~10x margin and
# queues far sooner than the 12 h "normal" maximum.
WALLTIME=00:30:00

CODE_DIR=/iopsstor/scratch/cscs/rroessner/code/ESMFold2_DuoDesign
LAUNCHER="$CODE_DIR/ESMFold2_Launch.sh"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Local subdirectory to cd into before sbatch, so relative paths in the job's
# own output (ESMFold2-DuoDesign_<jobid>_{out,err}.txt) land beside the results.
LOCAL_DIR=ESMFold2

if [ ! -x "$LAUNCHER" ]; then
    echo "ERROR: launcher not found or not executable: $LAUNCHER" >&2
    exit 1
fi

submit_one() {
    local target="$1" outdir="$2"
    shift 2
    local csv="$SCRIPT_DIR/${target}.csv"
    if [ ! -f "$csv" ]; then
        echo "ERROR: no design table for '$target' (expected $csv)" >&2
        return 1
    fi
    ( cd "$SCRIPT_DIR/$LOCAL_DIR" || exit 1
      set -x
      sbatch --time "$WALLTIME" --job-name "ESMFold2-$target" "$LAUNCHER" \
          "$csv" "./$outdir" --target "$target" --num-seeds "$NUM_SEEDS" \
          "$@" "${EXTRA_ARGS[@]}" )
}

mkdir -p "$SCRIPT_DIR/$LOCAL_DIR"

rc=0
for target in "${MSA_TARGETS[@]}"; do
    submit_one "$target" "$target" || rc=1
done
for target in "${NOMSA_TARGETS[@]}"; do
    submit_one "$target" "${target}_noMSA" --no-msa || rc=1
done

exit $rc
