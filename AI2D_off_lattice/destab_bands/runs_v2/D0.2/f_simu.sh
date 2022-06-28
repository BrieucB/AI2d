#!/bin/bash
#SBATCH --job-name=st_D0.2
#SBATCH -t 7-00:00:00
#SBATCH -n 16
#SBATCH --partition=multix
#SBATCH --nodelist=phoenix3,phoenix4
hostname
export OMP_NUM_THREADS=16

srun nice -n 19 /users/invites/benvegnen/Thesis/AI2D/AI2D_off_lattice/destab_bands/activeIsing

