#!/bin/bash
#SBATCH -o /home2/savila/HALOGEN/stderr.runs
##SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --export=NONE
#SBATCH --time=167:00:00
#SBATCH -A 16cores
##SBATCH --mem=600000
export OMP_NUM_THREADS=16
ulimit -c unlimited
time /home2/savila/libs/bin/mpiexec -np 16 2LPT-HALOGEN examples/FULL_HALOGEN/GOLIAT_2LPT.input examples/FULL_HALOGEN/GOLIAT_HALOGEN_run0.input
time /home2/savila/libs/bin/mpiexec -np 16 2LPT-HALOGEN examples/FULL_HALOGEN/GOLIAT_2LPT.input examples/FULL_HALOGEN/GOLIAT_HALOGEN_run1.input
time /home2/savila/libs/bin/mpiexec -np 16 2LPT-HALOGEN examples/FULL_HALOGEN/GOLIAT_2LPT.input examples/FULL_HALOGEN/GOLIAT_HALOGEN_run2.input
time /home2/savila/libs/bin/mpiexec -np 16 2LPT-HALOGEN examples/FULL_HALOGEN/GOLIAT_2LPT.input examples/FULL_HALOGEN/GOLIAT_HALOGEN_run3.input
time /home2/savila/libs/bin/mpiexec -np 16 2LPT-HALOGEN examples/FULL_HALOGEN/GOLIAT_2LPT.input examples/FULL_HALOGEN/GOLIAT_HALOGEN_run4.input
