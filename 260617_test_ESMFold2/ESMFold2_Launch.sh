#!/bin/bash -l

#SBATCH --partition gpu
#SBATCH --account dgfeller_epitope_pred
#SBATCH --job-name ESMFold2
#SBATCH --output %x_%j_out.txt
#SBATCH --error %x_%j_err.txt

#SBATCH --time 1-00:00:00
#SBATCH --mem 64G

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --gres gpu:A100:1
#SBATCH --gres-flags enforce-binding

#SBATCH --export NONE

# ESMFold2_Launch.sh
#
# SLURM job script to run ESMFold2 predictions via fold.py.
#
# Usage:
#   sbatch ESMFold2_Launch.sh <input.csv> <output_dir>
#   sbatch ESMFold2_Launch.sh -msa <input.csv> <output_dir> [msas_dir]
#
# To override SLURM defaults (e.g. time limit):
#   sbatch --time 2-00:00:00 ESMFold2_Launch.sh <input.csv> <output_dir>

USE_MSA=false
if [ "$1" == "-msa" ]; then
    USE_MSA=true
    shift
fi

if [ $# -lt 2 ]; then
    echo "Usage: sbatch ESMFold2_Launch.sh [-msa] <input.csv> <output_dir> [msas_dir]"
    exit 1
fi

INPUT_CSV=$(realpath "$1")
OUTPUT_DIR=$(realpath "$2")
MSAS_DIR=${3:+"$(realpath "$3")"}
SCRIPT_DIR=/scratch/rroessne/260617_test_ESMFold2

if [ "$USE_MSA" = true ]; then
    FOLD_SCRIPT="$SCRIPT_DIR/fold-msa.py"
    MSAS_DIR=${MSAS_DIR:-$(dirname "$INPUT_CSV")/msas}
else
    FOLD_SCRIPT="$SCRIPT_DIR/fold.py"
fi

echo -e "\n###########"
echo "[$(date +'%F %T')] Starting ESMFold2 predictions"
echo "  Input:  $INPUT_CSV"
echo "  Output: $OUTPUT_DIR"
echo "  Script: $FOLD_SCRIPT"
[ "$USE_MSA" = true ] && echo "  MSAs:   $MSAS_DIR"
echo "  Node:   $SLURM_JOB_NODELIST"
echo "  Job:    $SLURM_JOB_NAME ($SLURM_JOB_ID)"
echo -e "----------\n"

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

echo "Loading software stack..."
dcsrsoft use 20241118
echo "Loading miniforge3..."
module load miniforge3/24.11.3-2
echo "Running conda_init..."
conda_init
echo "Activating ESMFold2 env..."
conda activate ESMFold2
echo "Python: $(which python3)"

if [ "$USE_MSA" = true ]; then
    python3 -u "$FOLD_SCRIPT" "$INPUT_CSV" "$MSAS_DIR"
else
    python3 -u "$FOLD_SCRIPT" "$INPUT_CSV"
fi

echo -e "\n[$(date +'%F %T')] Finished ESMFold2 predictions."

exit 0
