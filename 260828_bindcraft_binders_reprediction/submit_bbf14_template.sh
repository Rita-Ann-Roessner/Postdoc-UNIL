#!/bin/bash -l

# submit_bbf14_template.sh -- launch BBF-14 DuoDesign folding jobs with no
# target MSA, for all 4 methods: ESMFold2 (plain --no-msa, no template --
# it's architecturally single-sequence, no MSA/template concept exists for
# it) and AF3/OpenDDE/Protenix (--no-msa plus the real 9HAC crystal
# structure manually supplied as a target template). Companion to
# submit_bbf14.sh (the plain no-msa-only, no-template comparison); run this
# for the "does a real structural template help/change folding without MSA"
# comparison.
#
# AF3 embeds the template mmCIF directly in its job JSON (see
# AF3_DuoDesign/build_input-msa.py's --template-cif). OpenDDE/Protenix use a
# different mechanism (proteinChain.templatesPath -> TemplateHitFeaturizer,
# see Protenix_DuoDesign/build_input-msa.py's --template-pdb-id docstring)
# that needed two extra pieces of infrastructure, both now in place:
#   - a `kalign` binary installed into the protenix/OpenDDE conda envs
#     (their template re-alignment step depends on it)
#   - a sitecustomize.py shim (DuoDesign_manual_template/, loaded via
#     PYTHONPATH by *_Launch.sh) that bypasses their template pipeline's
#     built-in duplicate/date prefilter, which would otherwise always
#     reject a hit this similar to the query -- see that file's docstring

# ─────────────────────────────────────────────────────────────────────────
# Which methods to actually submit. Comment out (or delete) any you don't
# want; order here doesn't matter. Valid names are the LOCAL_DIR keys below:
# ESMFold2, AF3, OpenDDE, Protenix. A name that isn't one of those is a hard
# error rather than a silent skip, so a typo can't quietly halve your run.
#METHODS=(
#    ESMFold2
#    AF3
#    OpenDDE
#    Protenix
#)
METHODS=(
    ESMFold2
)
# ─────────────────────────────────────────────────────────────────────────

TARGET=BBF-14

# Template source: the real 9HAC.cif (RCSB deposition of the BBF-14 /
# BBF-14_binder4 co-crystal structure).
#
# AF3 uses a hand-cleaned, chain-A-only extraction (see make_template_cif.py,
# which also carries over the release date AF3 insists on) -- AF3's own
# precomputed templates are single-chain-only mmCIF text, and it needs an
# explicit query<->template index mapping since it does no realignment of
# its own.
#
# Both index sets are 0-based *positional* indices, NOT residue numbers:
# BBF-14 residues 2-112 are query indices 1-111, and they map to the first
# 111 residues *present in the cif*, i.e. template indices 0-110 (the cif
# starts at label_seq_id 2 because Met1 is unresolved in the crystal -- that
# residue simply isn't in the file, so it does not occupy an index). Query
# index 0 (Met1) is therefore left unmapped, as are the cif's trailing
# GSG+His6 tag residues, which fall past index 110. Getting this wrong does
# NOT error -- AF3 folds against a frame-shifted template -- so verify with
# the featurisation snippet in AF3_DuoDesign/build_input-msa.py's docstring
# after changing any of these numbers. Verified 2026-09-20: this mapping
# reproduces the BBF-14 sequence at 111/112 positions (the 1 gap = Met1).
#
# OpenDDE/Protenix instead use the full original multi-chain 9HAC.cif
# (already placed at their expected mmcif-mirror locations -- see
# Protenix/params/mmcif/9hac.cif and cache/opendde/search_database/mmcif/9hac.cif)
# plus an explicit chain id: their own template pipeline re-aligns
# query<->template itself via kalign, so no manual index mapping is needed,
# and it already knows how to select one chain out of a multi-chain
# deposition via auth_chain_id filtering.
TEMPLATE_CIF="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/9HAC_chainA.cif"
TEMPLATE_QUERY_START=1
TEMPLATE_INDEX_START=0
TEMPLATE_LENGTH=111
TEMPLATE_PDB_ID=9HAC
TEMPLATE_CHAIN_ID=A

# 3 independent seeds, 1 sample per seed -- forwarded to every method's
# fold-msa.py (and, for AF3, also to build_input-msa.py's modelSeeds) via the
# shared --num-seeds flag supported by all four *_Launch.sh scripts.
NUM_SEEDS=3

DUODESIGN_ROOT=/work/FAC/FBM/LLB/dgfeller/epitope_pred/rroessne

# Local subdirectory (relative to this script) to cd into before sbatch, so
# ../BBF-14.csv and ./BBF-14_noMSA[_template] resolve correctly. Matches
# submit.sh's/submit_bbf14.sh's LOCAL_DIR mapping (Protenix_v2 is
# Protenix_Launch.sh's default model since no model_name is passed below).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
declare -A LOCAL_DIR=(
    [AF3]=AF3
    [ESMFold2]=ESMFold2
    [OpenDDE]=OpenDDE
    [Protenix]=Protenix_v2
)

if [ ${#METHODS[@]} -eq 0 ]; then
    echo "ERROR: METHODS is empty -- nothing to submit." >&2
    exit 1
fi
for m in "${METHODS[@]}"; do
    if [ -z "${LOCAL_DIR[$m]+set}" ]; then
        echo "ERROR: unknown method '$m' in METHODS." >&2
        echo "       Valid names: ${!LOCAL_DIR[*]}" >&2
        exit 1
    fi
done

should_run() {
    local m
    for m in "${METHODS[@]}"; do
        [ "$m" = "$1" ] && return 0
    done
    return 1
}

echo "Submitting ${TARGET}: ${METHODS[*]}"

# ESMFold2: plain --no-msa, no template. Same output dir convention
# (BBF-14_noMSA) as submit_bbf14.sh's other methods.
if should_run ESMFold2; then
(
    cd "$SCRIPT_DIR/${LOCAL_DIR[ESMFold2]}" || exit 1
    sbatch --time 01:30:00 "$DUODESIGN_ROOT/ESMFold2_DuoDesign/ESMFold2_Launch.sh" \
        "../${TARGET}.csv" "./${TARGET}_noMSA" --target "${TARGET}" --no-msa --num-seeds "${NUM_SEEDS}"
)
fi

# AF3: --no-msa + the manually supplied 9HAC template. Separate output dir
# (BBF-14_noMSA_template) from the plain no-msa run so both are directly
# comparable side by side.
if should_run AF3; then
(
    cd "$SCRIPT_DIR/${LOCAL_DIR[AF3]}" || exit 1
    sbatch --time 01:30:00 "$DUODESIGN_ROOT/AF3_DuoDesign/AF3_Launch.sh" \
        "../${TARGET}.csv" "./${TARGET}_noMSA_template" --target "${TARGET}" --no-msa --num-seeds "${NUM_SEEDS}" \
        --template-cif "${TEMPLATE_CIF}" \
        --template-query-start "${TEMPLATE_QUERY_START}" \
        --template-index-start "${TEMPLATE_INDEX_START}" \
        --template-length "${TEMPLATE_LENGTH}"
)
fi

# OpenDDE: --no-msa + the manually supplied 9HAC template.
if should_run OpenDDE; then
(
    cd "$SCRIPT_DIR/${LOCAL_DIR[OpenDDE]}" || exit 1
    sbatch --time 01:30:00 "$DUODESIGN_ROOT/OpenDDE_DuoDesign/OpenDDE_Launch.sh" \
        "../${TARGET}.csv" "./${TARGET}_noMSA_template" --target "${TARGET}" --no-msa --num-seeds "${NUM_SEEDS}" \
        --template-pdb-id "${TEMPLATE_PDB_ID}" --template-chain-id "${TEMPLATE_CHAIN_ID}"
)
fi

# Protenix: --no-msa + the manually supplied 9HAC template.
if should_run Protenix; then
(
    cd "$SCRIPT_DIR/${LOCAL_DIR[Protenix]}" || exit 1
    sbatch --time 01:30:00 "$DUODESIGN_ROOT/Protenix_DuoDesign/Protenix_Launch.sh" \
        "../${TARGET}.csv" "./${TARGET}_noMSA_template" --target "${TARGET}" --no-msa --num-seeds "${NUM_SEEDS}" \
        --template-pdb-id "${TEMPLATE_PDB_ID}" --template-chain-id "${TEMPLATE_CHAIN_ID}"
)
fi
