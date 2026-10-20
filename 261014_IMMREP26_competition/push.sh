#!/bin/bash -l
#
# push.sh -- copy the per-epitope TCR tables up to clariden.
#
# Sends <epitope>/<epitope>.csv from the local ESMFold2_predictions tree to the
# matching ESMFold2_predictions_clariden tree on the cluster. Remote per-epitope
# directories are created first in a single ssh call, since scp will not make
# them and fails per-file otherwise.
#
# Usage:  ./push.sh                  # every epitope in EPITOPES below
#         ./push.sh GILGFVFTL ...    # just these
#         DRY_RUN=1 ./push.sh        # print what would be copied

set -uo pipefail

EPITOPES=(CINGVCWTV GILGFVFTL GLCTLVAML LLWNGPMAV LTDEMIAQY NLVPMVATV TTDPSFLGRY YLQPRTFLL)

REMOTE=rroessner@clariden
LOCAL_ROOT=/Users/roessner/Documents/PostDoc/Data/261014_IMMREP26_competition/ESMFold2_predictions
REMOTE_ROOT=/iopsstor/scratch/cscs/rroessner/261014_IMMREP26_competition/ESMFold2_predictions_clariden

DRY_RUN="${DRY_RUN:-}"
(($# > 0)) && EPITOPES=("$@")

# Fail before touching the network if a table is missing locally.
missing=()
for ep in "${EPITOPES[@]}"; do
    [[ -f "$LOCAL_ROOT/$ep/$ep.csv" ]] || missing+=("$ep")
done
if ((${#missing[@]} > 0)); then
    echo "missing local tables: ${missing[*]}" >&2
    exit 1
fi

# One ssh round trip for all the directories rather than one per epitope.
mkdirs=""
for ep in "${EPITOPES[@]}"; do
    mkdirs+="mkdir -p '$REMOTE_ROOT/$ep'; "
done

if [[ -n "$DRY_RUN" ]]; then
    echo "ssh $REMOTE \"$mkdirs\""
else
    ssh "$REMOTE" "$mkdirs" || { echo "could not create remote directories" >&2; exit 1; }
fi

failed=()
for ep in "${EPITOPES[@]}"; do
    src="$LOCAL_ROOT/$ep/$ep.csv"
    dst="$REMOTE:$REMOTE_ROOT/$ep/$ep.csv"

    if [[ -n "$DRY_RUN" ]]; then
        echo "scp $src $dst"
        continue
    fi

    printf '%-12s %6s rows  -> %s\n' "$ep" "$(($(wc -l < "$src") - 1))" "$dst"
    scp -q "$src" "$dst" || failed+=("$ep")
done

if ((${#failed[@]} > 0)); then
    echo "FAILED: ${failed[*]}" >&2
    exit 1
fi
