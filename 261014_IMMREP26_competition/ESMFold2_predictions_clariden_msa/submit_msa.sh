#!/bin/bash -l
#
# submit_msa.sh -- the MSA arm of the IMMREP26 ESMFold2 runs.
#
# Same folds as ../ESMFold2_predictions_clariden, but with MSAs. Unlike
# templates (which ESMFold2 cannot take at all -- no released checkpoint
# carries the distogram-conditioning weights), MSAs are a NATIVE input: the
# model has an msa_encoder and fold-msa.py passes one per chain.
#
# THE THREE ARMS
#   control  ../ESMFold2_predictions_clariden                Fast,  no MSA
#   full     ../ESMFold2_predictions_clariden_full_checkpoint full,  no MSA
#   this     .                                               full, +MSA
#
# fold-msa.py hardcodes biohub/ESMFold2 (the full checkpoint) because Fast has
# msa_encoder.enabled=False and would silently ignore the MSAs. So this arm
# differs from the Fast control in TWO ways; the full arm is the control that
# isolates the MSA effect on its own.
#
# WHERE THE MSAs COME FROM
# No MSA search is run. build_input-msa.py ASSEMBLES a per-TCR MSA from
# pre-computed per-germline-gene AF3 MSAs under $PREMSA_ROOT, stitching the V
# and J gene alignments. Only ~150 gene segments exist and they are already
# computed, so this is seconds of CPU, not a search. Coverage of the genes used
# here is 152/153 -- only TRAJ58 is absent (27 TCRs), and those fall back to a
# thinner TCRA MSA rather than failing.
#
# Typical depths: MHC ~60k, B2M ~28k, TCRB ~11k, TCRA ~6k sequences.
#
# NOTHING RUNS ON THE LOGIN NODE
# Prep is a SLURM job (PrepMSA_Launch.sh, CPU only), and the fold job waits on
# it with --dependency=afterok, exactly like FitPrompt -> ESMFold2 in the
# template arm. Submitting returns immediately.
#
# DISK
# MHC.a3m (24.6 MB) and B2M.a3m (6.9 MB) are byte-identical for every complex
# of an epitope -- 31.5 of the 35.8 MB per complex. build_input-msa.py writes
# the bytes once and hardlinks the rest, so a tree is ~11 GB rather than ~94,
# and the page cache serves them across the node's 4 fold processes instead of
# re-reading 24 MB per fold. (dedup_msas.py does the same to an older tree
# that was built before that change.)
#
# COST (measured: 17.9 s/fold, vs 14.0 s full no-MSA and 6.8 s Fast)
#   LTDEMIAQY   390 ->  29 min     NLVPMVATV  1098 ->  82 min
#   GLCTLVAML   576 ->  43 min     TTDPSFLGRY 1161 ->  87 min
#   LLWNGPMAV   687 ->  51 min     YLQPRTFLL  1257 ->  94 min
#   CINGVCWTV   729 ->  54 min     GILGFVFTL  2010 -> 150 min
#   total 7908 folds ~ 39 GPU-hours
#
# Seeds match the other arms (base 1000, 3 seeds) so all three pair per
# (ID, seed).
#
# Usage:  ./submit_msa.sh                 # every epitope
#         ./submit_msa.sh GILGFVFTL ...   # just these
#         DRY_RUN=1 ./submit_msa.sh       # build MSAs + inputs, submit nothing
#         PREP_ONLY=1 ./submit_msa.sh     # build the MSAs, do not fold
#         REPREP=1 ./submit_msa.sh        # rebuild MSAs even if already present

set -uo pipefail

EPITOPES=(CINGVCWTV GILGFVFTL GLCTLVAML LLWNGPMAV LTDEMIAQY NLVPMVATV TTDPSFLGRY YLQPRTFLL)

NUM_SEEDS=3
BASE_SEED=1000

METHOD_ROOT=/iopsstor/scratch/cscs/rroessner/code
PREP_LAUNCHER="$METHOD_ROOT/ESMFold2/PrepMSA_Launch.sh"
LAUNCHER="$METHOD_ROOT/ESMFold2/ESMFold2_Launch.sh"

# build_input-msa.py reads all three from the environment, defaulting to the
# Curnagl /work paths. On Clariden they are in scratch.
export MIXTCRVIZ_DATA=/iopsstor/scratch/cscs/rroessner/MixTCRviz/data_raw/CDR123
export ALLELES_TAB=/iopsstor/scratch/cscs/rroessner/preMSA/Alleles_tab.txt
export PREMSA_ROOT=/iopsstor/scratch/cscs/rroessner/preMSA/AF3_models

TIME=12:00:00

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DRY_RUN="${DRY_RUN:-}"
PREP_ONLY="${PREP_ONLY:-}"
REPREP="${REPREP:-}"

[ $# -gt 0 ] && EPITOPES=("$@")

for epitope in "${EPITOPES[@]}"; do
    run_dir="$SCRIPT_DIR/$epitope"
    src="$run_dir/$epitope.csv"
    if [ ! -f "$src" ]; then
        echo "[skip] $epitope: no $epitope.csv"
        continue
    fi

    echo "=== $epitope"
    n_raw=$(( $(wc -l < "$src") - 1 ))

    # Is the MSA tree already built and complete?
    need_prep=1
    if [ -f "$run_dir/input.csv" ] && [ -d "$run_dir/msas" ]; then
        n_in=$(( $(wc -l < "$run_dir/input.csv") - 1 ))
        n_msa=$(ls "$run_dir/msas" 2>/dev/null | wc -l)
        if [ "$n_in" -gt 0 ] && [ "$n_in" -eq "$n_msa" ] && [ -z "$REPREP" ]; then
            need_prep=0
            echo "   prep     : cached ($n_msa complexes, $(du -sh "$run_dir/msas" 2>/dev/null | cut -f1)); REPREP=1 to rebuild"
        fi
    fi

    echo "   complexes: ~$n_raw x $NUM_SEEDS seeds = ~$(( n_raw * NUM_SEEDS )) folds" \
         "(~$(( n_raw * NUM_SEEDS * 18 / 4 / 60 )) min on 4 GPUs)"

    if [ -n "$DRY_RUN" ]; then
        [ $need_prep -eq 1 ] && echo "   [dry-run] would prep on a compute node, then fold" \
                             || echo "   [dry-run] would fold with the cached MSAs"
        continue
    fi

    dep=""
    if [ $need_prep -eq 1 ]; then
        jid=$(sbatch --parsable "$PREP_LAUNCHER" \
                "$src" "$run_dir/input.csv" "$run_dir/msas") || {
            echo "   [skip] prep submission failed"; continue; }
        echo "   prep job : $jid"
        dep="--dependency=afterok:$jid"
    fi

    if [ -n "$PREP_ONLY" ]; then
        echo "   (PREP_ONLY -- not submitting the fold)"
        continue
    fi

    # -msa routes the launcher to fold-msa.py and passes the msas dir through.
    fid=$( cd "$run_dir" && sbatch --parsable $dep --time "$TIME" "$LAUNCHER" -msa \
            "$run_dir/input.csv" "$run_dir" "$run_dir/msas" \
            --num-seeds "$NUM_SEEDS" --base-seed "$BASE_SEED" )
    if [ -n "$fid" ]; then
        echo "   fold job : $fid${dep:+  (waits on $jid)}"
    else
        echo "   [warn] fold submission failed"
    fi
done
