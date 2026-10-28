#!/bin/bash -l
#
# submit_templates.sh -- the TEMPLATE-PROMPTED arm of the IMMREP26 ESMFold2 runs.
#
# Same folds as ../ESMFold2_predictions_clariden (the unprompted control), but
# with the pMHC pinned to a reference structure by a fitted pair prompt. Two
# stages per epitope, chained with a SLURM dependency:
#
#   1. FitPrompt_Launch.sh   1 GPU    template -> prompts/<EPITOPE>.pt
#   2. ESMFold2_Launch.sh    4 GPUs   folds with --pair-prompt
#
# Both times below are TOTALS for the whole epitope, not per complex.
#
# Stage 1 is ~25 min REGARDLESS of how many TCRs the epitope has: the fit never
# sees a TCR. It runs on the pMHC alone (381 tokens, identical for every TCR of
# the epitope) for FIT_STEPS steps at ~5 s/step. That is the amortisation -- one
# 25 min fit covers every fold of that epitope. It is skipped when the prompt
# already exists, so re-running only refolds.
#
# Stage 2 scales with the epitope:  n_complexes x NUM_SEEDS x 6.8 s / 4 GPUs
#   LTDEMIAQY  130 ->  390 folds ~11 min     NLVPMVATV  366 -> 1098 ~31 min
#   GLCTLVAML  192 ->  576 folds ~16 min     TTDPSFLGRY 387 -> 1161 ~33 min
#   LLWNGPMAV  229 ->  687 folds ~19 min     YLQPRTFLL  419 -> 1257 ~36 min
#   CINGVCWTV  243 ->  729 folds ~21 min     GILGFVFTL  670 -> 2010 ~57 min
# (6.8 s/fold is the measured mean over all 7908 folds of the control arm.)
#
# ESMFold2 has no template input (its distogram-conditioning slot has no
# released weights, and PocketConditioning is a silent no-op), so the reference
# geometry goes in through the trunk's pair state instead. See
# code/ESMFold2/pair_prompt.py for why, and what that does and does not pin:
# it constrains CB positions, NOT side chains.
#
# Usage:  ./submit_templates.sh GLCTLVAML        # one epitope
#         ./submit_templates.sh                  # every epitope with a template
#         DRY_RUN=1 ./submit_templates.sh GLCTLVAML
#         ./submit_templates.sh GIL NLV ...             # prefixes also work
#         FIT_ONLY=1 ./submit_templates.sh GLCTLVAML    # fit, do not fold
#         REFIT=1 ./submit_templates.sh GLCTLVAML       # refit even if cached

set -uo pipefail

# ---------------------------------------------------------------------------
# Per-epitope reference structures.
#   <epitope>="<pdb under TEMPLATE_DIR>|<pipeline-chain=template-chain,...>"
# Add a line here when a new template lands. The chain map is checked at fit
# time: the template chain must be a substring of the pipeline sequence, so a
# wrong mapping fails immediately rather than fitting to the wrong thing.
# ---------------------------------------------------------------------------
declare -A TEMPLATE=(
  [CINGVCWTV]="template_pMHC_CINGVCWTV_3MRG.pdb|MHC=A,PEPTIDE=C"
  [GILGFVFTL]="template_pMHC_GILGFVFTL_1OGA.pdb|MHC=A,PEPTIDE=C"
  [GLCTLVAML]="template_pMHC_GLCTLVAML_3O4L.pdb|MHC=A,PEPTIDE=C"
  [LLWNGPMAV]="template_pMHC_LLWNGPMAV_5N6B.pdb|MHC=A,PEPTIDE=C"
  [NLVPMVATV]="template_pMHC_NLVPMVATV_3GSN.pdb|MHC=A,PEPTIDE=C"
  # Renamed 2026-10-21 from template_pMHC_NLVPMVATV_7N6E.pdb, which was
  # misnamed: 7N6E is the SARS-CoV-2 structure and its chain C is YLQPRTFLL.
  [YLQPRTFLL]="template_pMHC_YLQPRTFLL_7N6E.pdb|MHC=A,PEPTIDE=C"
)

# All six epitopes here are HLA-A*02:01 and every template's chain A is the
# same 180-residue alpha1-alpha2 groove; only the 9-mer differs. A new epitope
# therefore usually needs nothing but its own pMHC structure and one line above.

NUM_SEEDS=3
BASE_SEED=1000          # same seeds as the control arm, so the arms are paired

FIT_STEPS=300
FIT_SCOPE=peptide       # intra-peptide + peptide x groove only; the groove is
                        # already predicted to ~0.37 Ang unaided

METHOD_ROOT=/iopsstor/scratch/cscs/rroessner/code
BUILD_INPUT="$METHOD_ROOT/ESMFold2/build_input.py"
LAUNCHER="$METHOD_ROOT/ESMFold2/ESMFold2_Launch.sh"
FIT_LAUNCHER="$METHOD_ROOT/ESMFold2/FitPrompt_Launch.sh"
CHECK_TPL="$METHOD_ROOT/ESMFold2/template_peptide.py"
TEMPLATE_DIR=/iopsstor/scratch/cscs/rroessner/261014_IMMREP26_competition/templates

export MIXTCRVIZ_DATA=/iopsstor/scratch/cscs/rroessner/MixTCRviz/data_raw/CDR123
export ALLELES_TAB=/iopsstor/scratch/cscs/rroessner/preMSA/Alleles_tab.txt

TIME=12:00:00

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROMPT_DIR="$SCRIPT_DIR/prompts"
DRY_RUN="${DRY_RUN:-}"
FIT_ONLY="${FIT_ONLY:-}"
REFIT="${REFIT:-}"

if [ $# -gt 0 ]; then
    # Accept unambiguous prefixes (GIL -> GILGFVFTL) as well as full peptides.
    EPITOPES=()
    for arg in "$@"; do
        if [ -n "${TEMPLATE[$arg]:-}" ]; then
            EPITOPES+=("$arg"); continue
        fi
        hits=()
        for k in "${!TEMPLATE[@]}"; do
            case "$k" in "$arg"*) hits+=("$k");; esac
        done
        case ${#hits[@]} in
            1) EPITOPES+=("${hits[0]}") ;;
            0) echo "[skip] '$arg': no registered epitope starts with that" ;;
            *) echo "[skip] '$arg' is ambiguous: ${hits[*]}" ;;
        esac
    done
    [ ${#EPITOPES[@]} -eq 0 ] && exit 1
else
    EPITOPES=($(printf '%s\n' "${!TEMPLATE[@]}" | sort))
fi

mkdir -p "$PROMPT_DIR"

for epitope in "${EPITOPES[@]}"; do
    echo "=== $epitope"
    spec="${TEMPLATE[$epitope]:-}"
    if [ -z "$spec" ]; then
        echo "   [skip] no template registered -- add it to TEMPLATE above"
        continue
    fi
    run_dir="$SCRIPT_DIR/$epitope"
    src="$run_dir/$epitope.csv"
    if [ ! -f "$src" ]; then
        echo "   [skip] no $epitope.csv in $run_dir"
        continue
    fi

    tpl="$TEMPLATE_DIR/${spec%%|*}"
    chains="${spec##*|}"
    if [ ! -f "$tpl" ]; then
        echo "   [skip] template not found: $tpl"
        continue
    fi

    # Does the template really hold this epitope? Filenames lie (see the 7N6E
    # note above), and a template fitted to the wrong peptide would be a
    # plausible-looking prompt that pins the wrong thing. fit_pair_prompt.py
    # would reject it too, but failing here costs a second instead of a queue
    # wait plus a GPU allocation.
    pep_chain="${chains#*PEPTIDE=}"; pep_chain="${pep_chain%%,*}"
    tpl_pep=$(python3 "$CHECK_TPL" "$tpl" "$pep_chain" 2>/dev/null)
    want_pep=$(python3 -c "import csv,sys;print(next(iter(csv.DictReader(open(sys.argv[1]))))['peptide'])" "$src" 2>/dev/null)
    if [ -z "$tpl_pep" ] || [ "$tpl_pep" != "$want_pep" ]; then
        echo "   [skip] template peptide mismatch: $(basename "$tpl") chain"
        echo "          $pep_chain holds '$tpl_pep' but $epitope needs '$want_pep'"
        continue
    fi

    if ! python3 "$BUILD_INPUT" "$src" "$run_dir/input.csv"; then
        echo "   [skip] build_input.py failed"
        continue
    fi
    n=$(( $(wc -l < "$run_dir/input.csv") - 1 ))
    prompt="$PROMPT_DIR/$epitope.pt"

    echo "   template : $(basename "$tpl")  ($chains)"
    echo "   complexes: $n x $NUM_SEEDS seeds = $(( n * NUM_SEEDS )) folds"
    echo "   prompt   : $prompt"

    need_fit=1
    if [ -f "$prompt" ] && [ -z "$REFIT" ]; then
        need_fit=0
        echo "   fit      : cached (REFIT=1 to refit)"
    fi

    if [ -n "$DRY_RUN" ]; then
        [ $need_fit -eq 1 ] && echo "   [dry-run] would fit, then fold" \
                            || echo "   [dry-run] would fold with the cached prompt"
        continue
    fi

    dep=""
    if [ $need_fit -eq 1 ]; then
        jid=$(sbatch --parsable "$FIT_LAUNCHER" "$tpl" "$chains" \
                "$run_dir/input.csv" "$prompt" \
                --steps "$FIT_STEPS" --scope "$FIT_SCOPE") || {
            echo "   [skip] fit submission failed"; continue; }
        echo "   fit job  : $jid"
        dep="--dependency=afterok:$jid"
    fi

    if [ -n "$FIT_ONLY" ]; then
        echo "   (FIT_ONLY set -- not submitting the fold)"
        continue
    fi

    fid=$( cd "$run_dir" && sbatch --parsable $dep --time "$TIME" "$LAUNCHER" \
            "$run_dir/input.csv" "$run_dir" \
            --num-seeds "$NUM_SEEDS" --base-seed "$BASE_SEED" \
            --pair-prompt "$prompt" )
    if [ -n "$fid" ]; then
        echo "   fold job : $fid${dep:+  (waits on $jid)}"
    else
        echo "   [warn] fold submission failed"
    fi
done
