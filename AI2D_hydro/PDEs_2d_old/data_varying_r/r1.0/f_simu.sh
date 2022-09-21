#!/bin/bash
#SBATCH --job-name=pde1d
#SBATCH -t 13-00:00:00
#SBATCH -n 1
#SBATCH --nice=19
#SBATCH --partition=normal

hostname

srun nice -n 19 /users/invites/benvegnen/Thesis/AI1D/PDEs_1d/solve_PDEs_1d

