#!/bin/bash
#SBATCH --job-name=LES3
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH --account=fc_anemos
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=2400M
#SBATCH --time=71:59:59
#SBATCH --mail-type=all
#SBATCH --mail-user=adac@berkeley.edu

module purge
source ../../bashWRFV3

echo "=== Loaded Modules ==="
module list 2>&1
echo "======================"

mpiexec -n $SLURM_NTASKS ./ideal.exe
mpiexec -n $SLURM_NTASKS ./wrf.exe
