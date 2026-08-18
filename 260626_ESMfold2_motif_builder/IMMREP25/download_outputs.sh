#!/usr/bin/env bash
#
# Download a step's output.csv files (alpha + beta) from the cluster into the
# matching local atlas folders.
#
# Atlas layout has NO grid-point/condition level (unlike the benchmark):
#   Remote:  <REMOTE_BASE>/<epitope>/step<N>/<chain>/output.csv
#   Local:   <epitope>/step<N>/output_<chain>.csv
#
# Run from this atlas directory:  ./download_outputs.sh [step_number]   (default 1)

set -euo pipefail

REMOTE=rroessne@curnagl
REMOTE_BASE=/scratch/rroessne/260626_ESMfold2_motif_builder/IMMREP25

# Which step's outputs to fetch. Pass the step number as the first arg
# (e.g. ./download_outputs.sh 1); defaults to 1.
step=step2

# cd to the directory the script lives in (the local atlas root)
cd "$(dirname "$0")"

for epi in */; do
    epi=${epi%/}
    [[ -d "$epi/$step" ]] || continue
	
    for chain in alpha beta; do
        remote_path="$REMOTE_BASE/$epi/$step/$chain/output.csv"
        local_path="$epi/$step/output_${chain}.csv"

        echo ">> $chain: $epi"
        if scp "$REMOTE:$remote_path" "$local_path"; then
            echo "   ok -> $local_path"
        else
            echo "   FAILED: $remote_path" >&2
        fi

	remote_path="$REMOTE_BASE/$epi/$step/model_${chain}.csv"
        local_path="$epi/$step/input_${chain}.csv"

        echo ">> $chain: $epi"
        if scp "$REMOTE:$remote_path" "$local_path"; then
            echo "   ok -> $local_path"
        else
            echo "   FAILED: $remote_path" >&2
        fi

    done
done

echo "Done."
