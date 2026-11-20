#!/bin/bash
# Submit ESMFold2-Fast predictions for the listed CSVs, reduced cycles
# (num_loops=3, num_sampling_steps=50).
#
# These CSVs ALREADY contain the complete sequences (ID,TCRA,TCRB,MHC,B2M,PEPTIDE),
# so they go straight to fold.py -- no build_input.py step.
#
# Per CSV <name>.csv -> output dir ./<name>/ containing:
#   input.csv, output.csv, models/, metrics/, resources.csv, gpu_usage.csv

RUN_DIR=/iopsstor/scratch/cscs/rroessner/struc_preds_dummy
LAUNCHER=/iopsstor/scratch/cscs/rroessner/code/ESMFold2/ESMFold2_Launch.sh

NUM_LOOPS=3
NUM_SAMPLING_STEPS=50

CSVS=(
    LAU155_SLLMWITQC.csv
    LAU50_SLLMWITQC.csv
    motif_builder.csv
)

cd "$RUN_DIR" || exit 1

for csv in "${CSVS[@]}"; do
    if [ ! -f "$csv" ]; then
        echo "[skip] $csv not found"
        continue
    fi
    name="${csv%.csv}"
    echo "[submit] $name"
    sbatch "$LAUNCHER" "$RUN_DIR/$csv" "$RUN_DIR/$name" \
        --num-loops "$NUM_LOOPS" --num-sampling-steps "$NUM_SAMPLING_STEPS"
done
