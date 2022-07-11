#!/bin/bash
#SBATCH --job-name=tdoff_0.3
#SBATCH -t 13-00:00:00
#SBATCH -n 16
#SBATCH --partition=normalx
hostname
export OMP_NUM_THREADS=16

srun nice -n 19 /users/invites/benvegnen/Thesis/AI2D/AI2D_off_lattice/critical_D/activeIsing

