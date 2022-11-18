#!/bin/bash
#SBATCH --job-name=bdsD0.6
#SBATCH -t 7-00:00:00
#SBATCH --partition=multix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=8
hostname

export OMP_NUM_THREADS=8

srun nice -n 19 /users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/diag_rho_D_AIMrates/activeIsing

