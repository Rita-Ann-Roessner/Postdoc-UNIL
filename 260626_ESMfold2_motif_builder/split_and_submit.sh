#!/bin/bash
# =============================================================================
# Split a *_seqs.csv into batches of at most N TCRs and submit one ESMFold job
# per batch (so they run in parallel instead of one long serial job).
#
# Usage:  ./split_and_submit.sh <seqs.csv> <prefix> [N]
#   e.g.  ./split_and_submit.sh model_alpha_seqs.csv alpha 50
#
# Creates:   <prefix>_batches/batch_000.csv, batch_001.csv, ...
# Submits:   sbatch <LAUNCH> <batch.csv> ./<prefix>_batches/out_000  (per batch)
# Afterwards merge with:  ./merge_outputs.sh <prefix>
# =============================================================================
set -euo pipefail

LAUNCH=/lus/work/CT7/cad18049/roessner/ESMFold2/ESMFold2_Launch-Adastra.sh

infile=${1:?usage: split_and_submit.sh <seqs.csv> <prefix> [N]}
prefix=${2:?usage: split_and_submit.sh <seqs.csv> <prefix> [N]}
N=${3:-50}

[[ -f "$infile" ]] || { echo "no such file: $infile" >&2; exit 1; }

workdir="${prefix}_batches"
mkdir -p "$workdir"
hdr=$(head -1 "$infile")

# Split the data rows (header stripped) into N-row pieces: batch_000, batch_001...
tail -n +2 "$infile" | split -d -a 3 -l "$N" - "$workdir/batch_"

n=0
for f in "$workdir"/batch_[0-9]*; do
  [[ -e "$f" ]] || continue
  csv="$f.csv"
  { echo "$hdr"; cat "$f"; } > "$csv"   # re-add header to each chunk
  rm "$f"
  idx=$(basename "$f" | sed 's/batch_//')
  outdir="$workdir/out_$idx"
  sbatch "$LAUNCH" "$csv" "./$outdir"
  n=$((n + 1))
done

echo "submitted $n batch job(s) for $infile (<= $N TCRs each) -> $workdir/out_*"
