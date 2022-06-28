#!/bin/bash
#SBATCH --job-name=st_D0.7
#SBATCH -t 7-00:00:00
#SBATCH -n 16
#SBATCH --partition=multix
hostname
export OMP_NUM_THREADS=16

srun nice -n 19 /users/invites/benvegnen/Thesis/AI2D/AI2D_off_lattice/destab_bands/activeIsing

