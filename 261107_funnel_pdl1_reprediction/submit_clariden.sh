#!/bin/bash -l

# submit_clariden.sh -- repredictions of the PD-L1 funnel binders on
# Alps/Clariden. Every design in PDL1.csv is folded against BOTH targets by
# every method in METHODS: one sbatch job per (method, target).
#
#   hPDL1   human PD-L1, with its precomputed target MSA
#   mPDL1   mouse PD-L1, with its precomputed target MSA
#
# Same binder table every time, so the human-vs-mouse delta per design is the
# point -- a cross-species specificity read, not separate screens. Within a
# method, everything except the target is held fixed.
#
# Target MSAs live in code/ESMFold2_DuoDesign/preMSA/ and are shared by all
# four methods (each one's build_input-msa.py reads that same directory:
# ESMFold2/OpenDDE/Protenix take <Target>.a3m, AF3 takes the richer
# af3_jobs/<Target>_out/.../<target>_data.json). Binder chains never get an
# MSA -- they are de novo designs with no natural homologs.
#
# Each job gets one exclusive Clariden node (4x GH200); every *_Launch.sh
# splits that target's designs over the 4 GPUs and merges the results back.
#
# Run: ./submit_clariden.sh

# ─────────────────────────────────────────────────────────────────────────
# Which methods to submit. Valid names are the METHOD_DIR keys below.
# ESMFold2 is omitted by default because it has already been run to
# completion (ESMFold2_repredict/{h,m}PDL1, 2026-11-04); add it back to
# re-run it.
METHODS=(
    AF3
    OpenDDE
    Protenix-v2
)
# ─────────────────────────────────────────────────────────────────────────

# Targets to fold every design against, all WITH their precomputed MSA.
TARGETS=(hPDL1 mPDL1)

# 3 independent seeds, 1 sample per seed.
NUM_SEEDS=3

# Seeds are left unpinned: each per-GPU worker draws its own base seed, so a
# design gets NUM_SEEDS distinct seeds but the values differ between designs
# and between targets. The seeds used are recorded per fold in each run's
# output.csv. Add --base-seed <n> to a method's EXTRA_ARGS to align them.
#
# Per-method defaults are otherwise left exactly as each predictor ships them
# (ESMFold2: num_loops 20 / num_sampling_steps 200; AF3, OpenDDE and Protenix:
# cycle 10 / step 200 via --use_default_params), so nothing is passed here.
declare -A EXTRA_ARGS=(
    [ESMFold2]=""
    [AF3]=""
    [OpenDDE]=""
    [Protenix-v2]=""
)

# Wall time per job, now from measured per-fold cost on Clariden (4 designs x
# 1 seed, hPDL1, 2026-11-05) rather than the earlier guess of 4 h:
#   ESMFold2 ~15 s   AF3 ~24.5 s   OpenDDE ~20 s   Protenix ~17 s  per fold
# A full target is 30 designs x 3 seeds = 90 folds over 4 GPUs, i.e. ~23 folds
# per GPU: under 10 min of folding even for AF3, plus container start, model
# load and the MSA/JSON build. 1 h is therefore a ~4x margin for the three
# slower methods, and queues far sooner than a 4 h request.
declare -A WALLTIME=(
    [ESMFold2]=00:30:00
    [AF3]=01:00:00
    [OpenDDE]=01:00:00
    [Protenix-v2]=01:00:00
)

# Code dir holding each method's {method}_DuoDesign port, and the output dir
# prefix used under this directory. Protenix is "Protenix-v2" here because the
# model is protenix-v2 (Protenix_Launch.sh's own default), while its code dir
# is plain Protenix_DuoDesign.
declare -A METHOD_DIR=(
    [ESMFold2]=ESMFold2_DuoDesign
    [AF3]=AF3_DuoDesign
    [OpenDDE]=OpenDDE_DuoDesign
    [Protenix-v2]=Protenix_DuoDesign
)
declare -A LAUNCHER_NAME=(
    [ESMFold2]=ESMFold2_Launch.sh
    [AF3]=AF3_Launch.sh
    [OpenDDE]=OpenDDE_Launch.sh
    [Protenix-v2]=Protenix_Launch.sh
)

CODE_ROOT=/iopsstor/scratch/cscs/rroessner/code
PREMSA_DIR="$CODE_ROOT/ESMFold2_DuoDesign/preMSA"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BINDERS_CSV="$SCRIPT_DIR/PDL1.csv"

if [ ${#METHODS[@]} -eq 0 ]; then
    echo "ERROR: METHODS is empty -- nothing to submit." >&2
    exit 1
fi
for m in "${METHODS[@]}"; do
    if [ -z "${METHOD_DIR[$m]+set}" ]; then
        echo "ERROR: unknown method '$m' in METHODS." >&2
        echo "       Valid names: ${!METHOD_DIR[*]}" >&2
        exit 1
    fi
done
if [ ! -f "$BINDERS_CSV" ]; then
    echo "ERROR: no binder table at $BINDERS_CSV" >&2
    exit 1
fi

# Fail loudly if a target's MSA is missing, rather than letting a method fold
# it single-sequence and produce numbers that silently aren't comparable with
# the other target's. hPDL1 in particular was absent until it was ingested
# from preMSA/af3_jobs/hPDL1_out/ (precompute_target_msa.py ingest).
# ESMFold2/OpenDDE/Protenix need <Target>.a3m; AF3 needs the data.json.
for t in "${TARGETS[@]}"; do
    if [ ! -s "$PREMSA_DIR/${t}.a3m" ]; then
        echo "ERROR: no precomputed MSA for '$t' (expected $PREMSA_DIR/${t}.a3m)." >&2
        echo "       Produce it with:" >&2
        echo "         python3 $CODE_ROOT/ESMFold2_DuoDesign/precompute_target_msa.py make-job $t" >&2
        echo "       then run the AF3 data pipeline it prints and ingest the result." >&2
        exit 1
    fi
    for m in "${METHODS[@]}"; do
        [ "$m" = "AF3" ] || continue
        if ! compgen -G "$PREMSA_DIR/af3_jobs/${t}_out/*/*_data.json" > /dev/null; then
            echo "ERROR: AF3 needs $PREMSA_DIR/af3_jobs/${t}_out/*/*_data.json" >&2
            echo "       (the a3m is not enough -- AF3 embeds the MSA + templates itself)." >&2
            exit 1
        fi
    done
done

n_designs=$(( $(tail -n +2 "$BINDERS_CSV" | grep -c .) ))
echo "Binders: $BINDERS_CSV ($n_designs designs)"
echo "Methods: ${METHODS[*]}"
echo "Targets: ${TARGETS[*]}   seeds per design: $NUM_SEEDS (unpinned)"
echo

rc=0
for method in "${METHODS[@]}"; do
    launcher="$CODE_ROOT/${METHOD_DIR[$method]}/${LAUNCHER_NAME[$method]}"
    if [ ! -x "$launcher" ]; then
        echo "ERROR: launcher not found or not executable: $launcher" >&2
        rc=1
        continue
    fi
    work_dir="$SCRIPT_DIR/${method}_repredict"
    for target in "${TARGETS[@]}"; do
        mkdir -p "$work_dir/$target"
        # cd into the method's work dir so the job's own
        # <jobname>_<jobid>_{out,err}.txt land beside that method's results.
        ( cd "$work_dir" || exit 1
          set -x
          sbatch --time "${WALLTIME[$method]}" --job-name "${method}-${target}" "$launcher" \
              "$BINDERS_CSV" "./$target" --target "$target" \
              --num-seeds "$NUM_SEEDS" ${EXTRA_ARGS[$method]} ) || rc=1
    done
done

exit $rc
