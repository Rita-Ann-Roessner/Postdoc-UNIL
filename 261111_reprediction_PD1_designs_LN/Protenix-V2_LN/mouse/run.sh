#!/bin/bash
#SBATCH --job-name=nb_job
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=10:00:00
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=42G

mkdir -p logs

# Activate your environment
UMBRA_CONDA_ROOT_DEFAULT=""   # e.g. /work/lpdi/users/<you>/miniforge3
_conda_root="${UMBRA_CONDA_ROOT:-${UMBRA_CONDA_ROOT_DEFAULT}}"
if [[ -z "${_conda_root}" ]]; then
    command -v conda >/dev/null 2>&1 && _conda_root="$(conda info --base 2>/dev/null || true)"
    if [[ -z "${_conda_root}" ]]; then
        for _c in "${MAMBA_ROOT_PREFIX:-}" "${HOME}/miniforge3" "${HOME}/miniconda3" \
                  "${HOME}/mambaforge" "${HOME}/anaconda3" "/opt/conda"; do
            [[ -n "${_c}" && -x "${_c}/bin/activate" ]] && { _conda_root="${_c}"; break; }
        done
    fi
fi
readonly CONDA_ROOT="${_conda_root}"
readonly CONDA_ENV_NAME="${UMBRA_CONDA_ENV_NAME:-UmbraDesignTest2}"

if [[ -z "${CONDA_ROOT}" ]]; then
    echo "ERROR: could not locate a conda/mamba installation. Set UMBRA_CONDA_ROOT explicitly." >&2
    exit 1
fi

source "${CONDA_ROOT}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_NAME}"

# Sanity check — fail fast instead of silently running with the wrong/no env
python -c "import sys; print('Using interpreter:', sys.executable)"
echo "Active conda env: ${CONDA_DEFAULT_ENV}"

export UMBRA_REPO_ROOT="/scratch/nickel/UmbraOne_patched"   # or wherever UmbraCore/ actually lives

export UMBRA_REPO_ROOT="/scratch/nickel/UmbraOne_patched"   # or wherever UmbraCore/ actually lives
export PACKAGES_ROOT="/scratch/nickel/UmbraOne_patched"
export PROTENIX_ROOT_DIR="${PROTENIX_ROOT_DIR:-${PACKAGES_ROOT}/params/protenix}"
export PROTENIX_DATA_ROOT_DIR="${PROTENIX_DATA_ROOT_DIR:-${PROTENIX_ROOT_DIR}}"
jupyter nbconvert --to notebook --execute --inplace \
    /scratch/nickel/Reprediction_PD1_designs/Protenix-V2/mouse/protenix_reprediction.ipynb