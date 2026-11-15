#!/bin/bash -l

# submit_clariden.sh -- repredictions of the PD-1 designs on Alps/Clariden.
# Every design in PD1.csv is folded against BOTH targets by every method in
# METHODS: one sbatch job per (method, target).
#
#   hPD1    human PD-1, with its precomputed target MSA
#   mPD1    mouse PD-1, with its precomputed target MSA
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
# hPD1/mPD1 were not searched on Curnagl beforehand the way the PD-L1 targets
# were. They arrived as a prepared TargetPrep bundle (./Target/{hPD1,mPD1}/ --
# fasta, Protenix job spec, msa/, templates/), ingested into preMSA/ in both
# dialects with
#   code/ESMFold2_DuoDesign/ingest_targetprep_bundle.py ./Target hPD1 mPD1
# and then, for AF3 only, replaced by a real search run here:
#   sbatch code/ESMFold2_DuoDesign/AF3_DataPipeline_Launch.sh <target>
# (~8 min per target). ./Target itself is no longer in this directory; re-run
# the ingest from a fresh copy of the bundle if preMSA/ ever has to be rebuilt.
#
# So the two dialects in preMSA/ do NOT hold the same evidence:
#   <Target>.a3m               the bundle's ColabFold-style search. Deeper:
#                              4136 sequences for hPD1, 1607 for mPD1. Used by
#                              ESMFold2, OpenDDE and Protenix.
#   af3_jobs/<Target>_out/..   AF3's own data pipeline (max_template_date
#                              2021-09-30, as for PD-L1). Shallower unpaired
#                              MSA (325 kB vs 745 kB of rows for hPD1) but it
#                              also carries a paired MSA and four deposited PDB
#                              templates per target: 4 chains of 6ZQK at 47-83%
#                              coverage for hPD1, 4 chains of 7JG2 at 69-80%
#                              for mPD1. Used by AF3.
# The bundle's own single template -- a predicted self-template, pLDDT in its
# B-factor column -- is no longer used by anything.
#
# AF3 is also the only method here that can take a template at all: OpenDDE's
# template path wants a PDB ID already present in its own search_database,
# Protenix's needs the DuoDesign_manual_template shim (not set up on Alps), and
# ESMFold2 has no usable template channel. Those three fold the target from the
# MSA alone, exactly as they did for PD-L1.
#
# Each job gets one exclusive Clariden node (4x GH200); every *_Launch.sh
# splits that target's designs over the 4 GPUs and merges the results back.
#
# Run: ./submit_clariden.sh

# ─────────────────────────────────────────────────────────────────────────
# Which methods to submit. Valid names are the METHOD_DIR keys below.
# AF3 and OpenDDE already ran against both targets (jobs 3945262-3945265,
# 15 Nov 2026, 117/117 folds each); this pass adds the remaining two
# predictors, which read the same preMSA/<Target>.a3m OpenDDE did and so need
# no further preparation. Put all four names here to redo the whole screen.
METHODS=(
    ESMFold2
    Protenix-v2
)
# ─────────────────────────────────────────────────────────────────────────

# Targets to fold every design against, all WITH their precomputed MSA.
TARGETS=(hPD1 mPD1)

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

# Wall time per job, from measured per-fold cost on Clariden (hPDL1,
# 2026-11-05; the PD-1 targets are the same 115 aa, so the per-fold cost
# carries over):
#   ESMFold2 ~15 s   AF3 ~24.5 s   OpenDDE ~20 s   Protenix ~17 s  per fold
# A full target here is 39 designs x 3 seeds = 117 folds over 4 GPUs, i.e.
# ~30 folds per GPU: ~13 min of folding even for AF3, plus container start,
# model load and the MSA/JSON build. The whole PD-L1 screen (30 designs) ran
# in 3.5-9 min per job, so 1 h keeps a wide margin and queues far sooner than
# a larger request. AF3 additionally embeds the target MSA in every design
# JSON, which is ~0.9 MB per job for hPD1.
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
BINDERS_CSV="$SCRIPT_DIR/PD1.csv"

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
# the other target's -- build_input-msa.py only *warns* and falls back to
# single-sequence, which is why this is checked up front.
# ESMFold2/OpenDDE/Protenix need <Target>.a3m; AF3 needs the data.json.
for t in "${TARGETS[@]}"; do
    if [ ! -s "$PREMSA_DIR/${t}.a3m" ]; then
        echo "ERROR: no precomputed MSA for '$t' (expected $PREMSA_DIR/${t}.a3m)." >&2
        echo "       For a target prepared elsewhere (TargetPrep bundle -- note that" >&2
        echo "       ./Target was removed after the ingest, so copy it back first), run:" >&2
        echo "         python3 $CODE_ROOT/ESMFold2_DuoDesign/ingest_targetprep_bundle.py \\" >&2
        echo "             $SCRIPT_DIR/Target $t" >&2
        echo "       To search it here instead, use precompute_target_msa.py make-job $t," >&2
        echo "       run the AF3 data pipeline it prints, then ingest the result." >&2
        exit 1
    fi
    for m in "${METHODS[@]}"; do
        [ "$m" = "AF3" ] || continue
        if ! compgen -G "$PREMSA_DIR/af3_jobs/${t}_out/*/*_data.json" > /dev/null; then
            echo "ERROR: AF3 needs $PREMSA_DIR/af3_jobs/${t}_out/*/*_data.json" >&2
            echo "       (the a3m is not enough -- AF3 embeds the MSA + templates itself);" >&2
            echo "       ingest_targetprep_bundle.py writes both from one bundle." >&2
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
