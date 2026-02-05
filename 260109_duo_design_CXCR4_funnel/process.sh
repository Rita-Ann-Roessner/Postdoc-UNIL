#!/bin/bash
#SBATCH --job-name=MD-BindCraft
#SBATCH --nodes=1
################
#SBATCH --ntasks-per-node=1         # (number of mpi tasks) / 4 MPI par replica  
#SBATCH --cpus-per-task=8          # (Number of OpenMP threads per MPI task) / nombre de coeurs à réserver 
#SBATCH --gres=gpu:1                # Number of GPUs per node à réserver
#SBATCH --hint=nomultithread               
#SBATCH -C MI250                
#SBATCH -A c1615115
#SBATCH --exclusive # temporary until module error is fixed
######
#SBATCH --time=00:10:00
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

#-----------------

pwd=$(pwd)
parent_dir=$(dirname "$pwd")
grandparent_dir=$(dirname $(dirname "$pwd"))

echo q | srun -K1 --mpi=cray_shasta -N 1 -n 1 -c "$threads" -m block:block --cpu-bind cores \
gmx_mpi make_ndx -f npt.gro

echo System | srun -K1 --mpi=cray_shasta -N 1 -n 1 -c "$threads" -m block:block --cpu-bind cores \
gmx_mpi trjconv -f production.xtc -s npt.gro -n index.ndx -pbc nojump -o tmp.xtc -dt 1000

echo Protein Protein | srun -K1 --mpi=cray_shasta -N 1 -n 1 -c "$threads" -m block:block --cpu-bind cores \
gmx_mpi trjconv -f tmp.xtc -s production.tpr -n index.ndx -fit rot+trans -o reimage.xtc

echo Protein | srun -K1 --mpi=cray_shasta -N 1 -n 1 -c "$threads" -m block:block --cpu-bind cores \
gmx_mpi editconf -f npt.gro -n index.ndx -o reimage.gro 

rm \#*\#

