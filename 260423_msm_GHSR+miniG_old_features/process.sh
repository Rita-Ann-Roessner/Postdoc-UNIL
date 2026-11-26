#!/bin/bash
#SBATCH --job-name=MD-MSM
#SBATCH --nodes=1
################
#SBATCH --ntasks-per-node=1        # (number of mpi tasks) / 4 MPI par replica  
#SBATCH --cpus-per-task=10          # (Number of OpenMP threads per MPI task) / nombre de coeurs à réserver 
#SBATCH --gres=gpu:1                # Number of GPUs per node à réserver
#SBATCH --hint=nomultithread
#SBATCH --exclusive               
#SBATCH -C MI250                
#SBATCH -A cad18049
######
#SBATCH --time=2:00:00
#SBATCH -o job.out
#SBATCH -e job.err
######
module purge
module load archive CCE-GPU-3.1.0
module load gromacs/2023_amd-omp-plumed-rocm-mpi-gpu

module list

nodes="${SLURM_JOB_NUM_NODES}"
threads="${SLURM_CPUS_PER_TASK}"
ranks="${SLURM_NTASKS}"

# OpenMP thread binding
export OMP_NUM_THREADS="$threads"
export OMP_PLACES=cores
export OMP_PROC_BIND=close

# NIC POLICY
export MPICH_OFI_NIC_POLICY=GPU

# Stack limits
ulimit -c unlimited
ulimit -s unlimited      


# gromacs settings
export GMX_DISABLE_GPU_TIMING=1
export GMX_GPU_DD_COMMS=1
export GMX_GPU_PME_PP_COMMS=1
export GMX_NO_QUOTES=1
export GMX_ENABLE_DIRECT_GPU_COMM=1

-----------------

pwd=$(pwd)
parent_dir=$(dirname "$pwd")
grandparent_dir=$(dirname "$parent_dir")
INPUT="${grandparent_dir}/MD_input/"
export JOB="production"

for dir in {0..3}
do
    cd ${pwd}/${dir}
    echo SOLU System | srun -K1 --mpi=cray_shasta -N 1 -n 1 -c "$threads" -m block:block --cpu-bind cores gmx_mpi trjconv -f init.pdb -s ${JOB}.tpr -n ${INPUT}/index.ndx -ur compact -pbc mol -center -o reimage.gro
    echo SOLU System | srun -K1 --mpi=cray_shasta -N 1 -n 1 -c "$threads" -m block:block --cpu-bind cores gmx_mpi trjconv -f ${JOB}.xtc -s ${JOB}.tpr -n ${INPUT}/index.ndx -ur compact -pbc mol -center -o reimage.xtc

    echo SOLU SOLU | srun -K1 --mpi=cray_shasta -N 1 -n 1 -c "$threads" -m block:block --cpu-bind cores gmx_mpi trjconv -f init.pdb -s ${JOB}.tpr -n ${INPUT}/index.ndx -ur compact -pbc mol -center -o protein.gro
    echo SOLU SOLU | srun -K1 --mpi=cray_shasta -N 1 -n 1 -c "$threads" -m block:block --cpu-bind cores gmx_mpi trjconv -f ${JOB}.xtc -s ${JOB}.tpr -n ${INPUT}/index.ndx -ur compact -pbc mol -center -o protein.xtc

    cd ${pwd}
done

