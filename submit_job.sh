#!/bin/bash
#SBATCH --account=ucb598_asc1
#SBATCH --job-name=WRF-MITgcm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=3840M
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --time=24:00:00
#SBATCH --output=WRF-MITgcm.%j.log

source /home/yuli3660/projects/library-env-vars.sh

cd /projects/yuli3660/test_coarse/SKRIPS/runCase

echo "Executing mpirun -np ... ./esmf_application"
mpirun -np 4 ./esmf_application
