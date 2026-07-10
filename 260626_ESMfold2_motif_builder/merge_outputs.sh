#!/bin/bash
# =============================================================================
# Merge the per-batch ESMFold outputs produced by split_and_submit.sh back into
# ONE directory that reproduces the original single-job layout:
#     <prefix>/output.csv          (all batch output.csv, one header)
#     <prefix>/models/tcr*.cif     (all structures, gathered)
#     <prefix>/metrics/tcr*.csv    (all per-model metrics, gathered)
# Run after all batch jobs finish.
#
# Usage:  ./merge_outputs.sh <prefix>
#   e.g.  ./merge_outputs.sh alpha   ->  alpha/{output.csv, models, metrics}
#
# Files are named by TCR id (tcr<id>.cif / tcr<id>.csv), which is globally unique
# across batches, so gathering is collision-free. A collision is still checked
# for and reported (it would mean index-named files, not id-named).
#
# Afterwards, place <prefix>/output.csv into the step dir your R pipeline reads
# (as output_<chain>.csv), exactly as with the old single-job workflow.
# =============================================================================
set -euo pipefail

prefix=${1:?usage: merge_outputs.sh <prefix>}
workdir="${prefix}_batches"
dest="$prefix"

[[ -d "$workdir" ]] || { echo "no batch dir: $workdir" >&2; exit 1; }

shopt -s nullglob
dirs=("$workdir"/out_*)
[[ ${#dirs[@]} -gt 0 ]] || { echo "no out_* dirs in $workdir" >&2; exit 1; }

mkdir -p "$dest/models" "$dest/metrics"

missing=0; collide=0; first=1
for d in "${dirs[@]}"; do
  o="$d/output.csv"
  if [[ ! -f "$o" ]]; then
    echo "  MISSING: $o (job not finished / failed?)" >&2
    missing=$((missing + 1)); continue
  fi

  # merged scores table (keep a single header)
  if [[ $first -eq 1 ]]; then cp "$o" "$dest/output.csv"; first=0
  else tail -n +2 "$o" >> "$dest/output.csv"; fi

  # gather structures + per-model metrics (copy, so re-runnable; collision-checked)
  for sub in models metrics; do
    [[ -d "$d/$sub" ]] || continue
    for fpath in "$d/$sub"/*; do
      [[ -e "$fpath" ]] || continue
      fname=$(basename "$fpath")
      if [[ -e "$dest/$sub/$fname" ]]; then
        echo "  COLLISION: $sub/$fname already present (index-named, not id-named?)" >&2
        collide=$((collide + 1))
      else
        cp "$fpath" "$dest/$sub/$fname"
      fi
    done
  done
done

[[ $missing -gt 0 ]] && \
  echo "WARNING: $missing/${#dirs[@]} batch(es) missing output.csv — INCOMPLETE" >&2
[[ $collide -gt 0 ]] && \
  echo "WARNING: $collide filename collision(s) — models/metrics NOT fully gathered" >&2

rows=$(( $(wc -l < "$dest/output.csv") - 1 ))
ncif=$(find "$dest/models"  -maxdepth 1 -name '*.cif' | wc -l)
nmet=$(find "$dest/metrics" -maxdepth 1 -name '*.csv' | wc -l)
echo "merged $(( ${#dirs[@]} - missing ))/${#dirs[@]} batches -> $dest/  (output.csv: $rows TCRs, $ncif cif, $nmet metrics)"
