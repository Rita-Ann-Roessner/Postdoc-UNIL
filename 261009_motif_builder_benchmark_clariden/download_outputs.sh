#!/usr/bin/env bash
#
# Download a step's CSV tables (alpha + beta) from the cluster into the
# matching local benchmark folders.
#
# Benchmark layout has a METHOD level on top (unlike the atlas), and on the
# cluster the per-chain results sit in their own subdirectory. Those are
# flattened on the way down, to the same output_<chain>.csv / input_<chain>.csv
# convention the atlas folders use:
#
#   Remote:  <REMOTE_BASE>/<method>/<epitope>/<step>/*.txt
#            <REMOTE_BASE>/<method>/<epitope>/<step>/model_<chain>.csv
#            <REMOTE_BASE>/<method>/<epitope>/<step>/model_<chain>_seqs.csv
#            <REMOTE_BASE>/<method>/<epitope>/<step>/<chain>/input.csv
#            <REMOTE_BASE>/<method>/<epitope>/<step>/<chain>/output.csv
#            <REMOTE_BASE>/<method>/<epitope>/<step>/<chain>/gpu_usage.csv
#            <REMOTE_BASE>/<method>/<epitope>/<step>/<chain>/resources.csv
#   Local:   <method>/<epitope>/<step>/*.txt
#            <method>/<epitope>/<step>/model_<chain>.csv
#            <method>/<epitope>/<step>/model_<chain>_seqs.csv
#            <method>/<epitope>/<step>/input_<chain>.csv
#            <method>/<epitope>/<step>/output_<chain>.csv
#            <method>/<epitope>/<step>/gpu_usage_<chain>.csv
#            <method>/<epitope>/<step>/resources_<chain>.csv
#
# No alpha/ or beta/ subdirectories locally, and nothing heavy comes down --
# the models/, metrics/ and *_batches/ trees stay on the cluster. A full run is
# about 18 MB.
#
# The per-chain input.csv is the sequence-level table that was actually folded
# (ID,TCRA,TCRB,MHC,B2M,PEPTIDE); model_<chain>.csv is the gene-level table it
# was built from. Both come down because the R pipeline keys on columns of each.
#
# Run from the directory you want the tree in:
#   ./download_outputs.sh                 # step0, all methods
#   ./download_outputs.sh step2           # a different step
#   ./download_outputs.sh step0 AF3       # one method
#   ./download_outputs.sh 2 AF3 Protenix  # bare step number also works
#   DRY_RUN=1 ./download_outputs.sh       # print what would be fetched, fetch nothing

set -uo pipefail

# Both overridable from the environment, e.g. if your ssh config names the
# cluster differently:  REMOTE=curnagl ./download_outputs.sh
REMOTE="${REMOTE:-rroessner@clariden}"
REMOTE_BASE="${REMOTE_BASE:-/iopsstor/scratch/cscs/rroessner/261009_motif_builder_benchmark_clariden}"

# Which step to fetch, as an optional first argument -- either "step2" or a
# bare "2"; defaults to step0. Anything else is taken as a method name rather
# than a step, so "./download_outputs.sh AF3" does what it looks like.
step=step0
if [[ $# -gt 0 ]]; then
    case "$1" in
        step[0-9]*)      step="$1";      shift ;;
        [0-9]|[0-9][0-9]) step="step$1"; shift ;;
    esac
fi

# Method directories to pull. Override by passing names after the step.
METHODS=(esmfold2)
[[ $# -gt 0 ]] && METHODS=("$@")

CHAINS=(alpha beta)

# Local folder name for a remote method dir, when the two differ. A case
# statement rather than an associative array, so macOS's bash 3.2 runs it.
local_name() {
    case "$1" in
        esmfold2) echo ESMFold2 ;;
        *)        echo "$1" ;;
    esac
}

# Epitope directories are normally discovered on the cluster. Set EPITOPES to a
# space-separated list to skip that lookup (useful offline, or if ssh is slow):
#   EPITOPES="A0201_LLWNGPMAV" ./download_outputs.sh
EPITOPES="${EPITOPES:-}"

# cd to the directory the script lives in (the local benchmark root)
cd "$(dirname "$0")"

DRY_RUN="${DRY_RUN:-}"

ok=0
failed=0

# Share one SSH connection across every ssh/scp below, so a password-protected
# key is asked for once instead of ~30 times. Harmless if it fails: SSH_OPTS
# stays empty and each transfer opens its own connection as usual.
SSH_OPTS=()
CTL="${TMPDIR:-/tmp}/.download_outputs_ctl_$$"
if ssh -o ControlMaster=yes -o ControlPersist=120 -o ControlPath="$CTL" \
       -fN "$REMOTE" 2>/dev/null; then
    SSH_OPTS=(-o ControlPath="$CTL")
    trap 'ssh -O exit -o ControlPath="$CTL" "$REMOTE" 2>/dev/null || true' EXIT
fi

# fetch <remote-relative-path> <local-path>
fetch() {
    local rel="$1" dest="$2"
    if [[ -n "$DRY_RUN" ]]; then
        echo "   [dry-run] $rel -> $dest"
        return 0
    fi
    mkdir -p "$(dirname "$dest")"
    if scp -q "${SSH_OPTS[@]}" "$REMOTE:$REMOTE_BASE/$rel" "$dest"; then
        echo "   ok -> $dest"
        ok=$((ok + 1))
    else
        echo "   FAILED: $rel" >&2
        failed=$((failed + 1))
    fi
}

for method in "${METHODS[@]}"; do
    # Discovered on the cluster, so this works on a fresh checkout where
    # nothing exists locally yet. A connection failure is reported separately
    # from an existing-but-empty method dir, so a broken ssh is never silently
    # read as "this method has no epitopes".
    ldir=$(local_name "$method")
    if [[ -n "$EPITOPES" ]]; then
        epitopes="$EPITOPES"
    else
        epitopes=$(ssh "${SSH_OPTS[@]}" "$REMOTE" \
            "cd '$REMOTE_BASE/$method' 2>/dev/null && ls -d */ 2>/dev/null" | tr -d '/')
        rc=$?
        if [[ $rc -eq 255 ]]; then
            echo "[error] $method: cannot reach $REMOTE (ssh failed). Check your" >&2
            echo "        ssh key / host alias, or set EPITOPES= to skip the lookup." >&2
            failed=$((failed + 1))
            continue
        fi
    fi
    if [[ -z "$epitopes" ]]; then
        echo "[skip] $method: no epitope dirs under $REMOTE_BASE/$method" >&2
        continue
    fi

    for epi in $epitopes; do
        echo ">> $method/$epi"

        # Whatever *.txt files sit at the step level; names vary, so they are
        # listed on the cluster first and then fetched one by one.
        txts=$(ssh "${SSH_OPTS[@]}" "$REMOTE" \
            "cd '$REMOTE_BASE/$method/$epi/$step' 2>/dev/null && ls *.txt 2>/dev/null")
        for txt in $txts; do
            fetch "$method/$epi/$step/$txt" "$ldir/$epi/$step/$txt"
        done

        for chain in "${CHAINS[@]}"; do
            echo ">> $chain: $method/$epi/$step"

            # Gene-level and sequence-level model tables keep their names.
            fetch "$method/$epi/$step/model_${chain}.csv" \
                  "$ldir/$epi/$step/model_${chain}.csv"
            fetch "$method/$epi/$step/model_${chain}_seqs.csv" \
                  "$ldir/$epi/$step/model_${chain}_seqs.csv"

            # The two per-chain tables get flattened out of <chain>/ and
            # renamed, so alpha and beta can share one step directory.
            fetch "$method/$epi/$step/$chain/input.csv" \
                  "$ldir/$epi/$step/input_${chain}.csv"
            fetch "$method/$epi/$step/$chain/output.csv" \
                  "$ldir/$epi/$step/output_${chain}.csv"
            fetch "$method/$epi/$step/$chain/gpu_usage.csv" \
                  "$ldir/$epi/$step/gpu_usage_${chain}.csv"
            fetch "$method/$epi/$step/$chain/resources.csv" \
                  "$ldir/$epi/$step/resources_${chain}.csv"
        done
    done
done

echo
if [[ -n "$DRY_RUN" ]]; then
    echo "DRY RUN: nothing was written."
else
    echo "Done. $ok file(s) fetched, $failed failed."
fi
