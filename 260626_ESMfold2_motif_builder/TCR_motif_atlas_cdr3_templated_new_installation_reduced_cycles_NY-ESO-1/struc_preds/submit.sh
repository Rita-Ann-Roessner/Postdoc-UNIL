#!/bin/bash
# Submit ESMFold2-Fast predictions for every TCR-table CSV in this directory,
# with reduced cycles (num_loops=3, num_sampling_steps=50).
#
# Each CSV is a TCR table (the step before sequences) with columns:
#   id, TRAV, cdr3_TRA, TRBV, cdr3_TRB, TRAJ, TRBJ, model, peptide, MHC_allele_a
# build_input.py turns that into fold.py's sequence input (ID,TCRA,TCRB,MHC,B2M,PEPTIDE);
# species is inferred from MHC_allele_a, the extra `model` column is ignored.
#
# Per CSV <name>.csv -> output dir ./<name>/ containing:
#   input.csv, output.csv, models/, metrics/, resources.csv, gpu_usage.csv

RUN_DIR=/iopsstor/scratch/cscs/rroessner/struc_preds
CODE=/iopsstor/scratch/cscs/rroessner/code

NUM_LOOPS=3
NUM_SAMPLING_STEPS=50

# Reference tables for build_input.py (live on Clariden scratch).
export MIXTCRVIZ_DATA=/iopsstor/scratch/cscs/rroessner/MixTCRviz/data_raw/CDR123
export ALLELES_TAB=/iopsstor/scratch/cscs/rroessner/preMSA/Alleles_tab.txt

cd "$RUN_DIR" || exit 1

for csv in *.csv; do
    [ -e "$csv" ] || { echo "[skip] no CSVs found"; break; }
    name="${csv%.csv}"
    echo "[submit] $name"
    mkdir -p "$RUN_DIR/$name"
    python3 "$CODE/ESMFold2/build_input.py" "$csv" "$RUN_DIR/$name/input.csv" || {
        echo "ERROR: build_input.py failed for $csv" >&2; continue; }
    sbatch "$CODE/ESMFold2/ESMFold2_Launch.sh" "$RUN_DIR/$name/input.csv" "$RUN_DIR/$name" \
        --num-loops "$NUM_LOOPS" --num-sampling-steps "$NUM_SAMPLING_STEPS"
done
