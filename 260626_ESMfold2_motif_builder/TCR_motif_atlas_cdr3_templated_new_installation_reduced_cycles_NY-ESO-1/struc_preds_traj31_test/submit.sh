#!/bin/bash
# ESMFold2-Fast predictions for the TRAV17/TRAJ31 test panel, reduced cycles
# (num_loops=3, num_sampling_steps=50) -- MUST match the step-0 settings, or
# the scores are not comparable to the panel this test is interpreted against.

RUN_DIR=/iopsstor/scratch/cscs/rroessner/struc_preds_traj31_test
LAUNCHER=/iopsstor/scratch/cscs/rroessner/code/ESMFold2/ESMFold2_Launch.sh

NUM_LOOPS=3
NUM_SAMPLING_STEPS=50

CSVS=(
    traj31_test.csv
)

cd "$RUN_DIR" || exit 1

for csv in "${CSVS[@]}"; do
    if [ ! -f "$csv" ]; then
        echo "[skip] $csv not found"
        continue
    fi
    n="${csv%.csv}"
    echo "[submit] $n"
    sbatch "$LAUNCHER" "$RUN_DIR/$csv" "$RUN_DIR/$n" \
        --num-loops "$NUM_LOOPS" --num-sampling-steps "$NUM_SAMPLING_STEPS"
done
