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
##SBATCH --exclusive # temporary until module error is fixed
######
#SBATCH --time=24:00:00
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
grandparent_dir=$(dirname "$parent_dir")

export SUBMIT="job.sh"
START=$(date +"%s")

if [ -e production.tpr ];then
    env OMP_NUM_THREADS=$threads srun -K1 --mpi=cray_shasta -N "$nodes" -n "$ranks" -c "$threads" -m block:block --cpu-bind cores \
    gmx_mpi mdrun -deffnm production -cpi production.cpt -ntomp $threads -plumed plumed.dat -maxh 23.9
else
    # nvt
    srun -K1 --mpi=cray_shasta -N 1 -n 1 -c "$threads" -m block:block --cpu-bind cores \
    gmx_mpi grompp -f ${grandparent_dir}/nvt.mdp -o nvt.tpr -c ${parent_dir}/minimization.gro -r ${parent_dir}/minimization.gro -p ${parent_dir}/topol.top -n ${parent_dir}/index.ndx -maxwarn 1

    env OMP_NUM_THREADS=$threads srun -K1 --mpi=cray_shasta -N "$nodes" -n "$ranks" -c "$threads" -m block:block --cpu-bind cores \
    gmx_mpi mdrun -deffnm nvt -ntomp $threads

    # npt
    srun -K1 --mpi=cray_shasta -N 1 -n 1 -c "$threads" -m block:block --cpu-bind cores \
    gmx_mpi grompp -f ${grandparent_dir}/npt.mdp -o npt.tpr -c nvt.gro -r nvt.gro -p ${parent_dir}/topol.top -n ${parent_dir}/index.ndx -maxwarn 2

    env OMP_NUM_THREADS=$threads srun -K1 --mpi=cray_shasta -N "$nodes" -n "$ranks" -c "$threads" -m block:block --cpu-bind cores \
    gmx_mpi mdrun -deffnm npt -ntomp $threads

    # production
    srun -K1 --mpi=cray_shasta -N 1 -n 1 -c "$threads" -m block:block --cpu-bind cores \
    gmx_mpi grompp -f ${grandparent_dir}/production.mdp -o production.tpr -c npt.gro -p ${parent_dir}/topol.top -n ${parent_dir}/index.ndx -maxwarn 2

    env OMP_NUM_THREADS=$threads srun -K1 --mpi=cray_shasta -N "$nodes" -n "$ranks" -c "$threads" -m block:block --cpu-bind cores \
    gmx_mpi mdrun -deffnm production -cpi production.cpt -ntomp $threads -plumed plumed.dat -maxh 20
fi


END=$(date +"%s")
echo "$(((END-START))) seconds ran"
echo "$(((END-START)/3600)) full hours ran"

num1=$(bc <<< "scale=4;($END - $START) / 3600")
num2=23

if (( $(echo "$num1 < $num2" |bc -l) )); then
  echo "last cycle was just $(((END-START)/60))min long and therefore finito"
  exit 3
else
  echo "cycle resubmitting"
  sbatch ${SUBMIT}
  exit 2
fi

