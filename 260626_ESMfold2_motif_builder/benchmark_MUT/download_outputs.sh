#!/usr/bin/env bash
#
# Download step1 output.csv files (alpha + beta) from the cluster into the
# matching local benchmark folders.
#
# Remote layout:  <REMOTE_BASE>/<cond>/<epitope>/step1/<chain>/output.csv
# Local  layout:  <cond>/<epitope>/step1/output_<chain>.csv
#
# Run from the benchmark/ directory:  ./download_outputs.sh

set -euo pipefail

REMOTE=rroessne@curnagl
REMOTE_BASE=/scratch/rroessne/260626_ESMfold2_motif_builder/benchmark_MUT

# Which step's outputs to fetch. 
step=step2

# cd to the directory the script lives in (the local benchmark/ root)
cd "$(dirname "$0")"

for cond in */; do
    cond=${cond%/}
    [[ -d "$cond" ]] || continue

    for epi in "$cond"/*/; do
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
        done
    done
done

echo "Done."
