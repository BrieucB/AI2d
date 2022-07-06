#!/bin/bash
#SBATCH --job-name=1n_ly
#SBATCH -t 7-00:00:00
#SBATCH -n 8
#SBATCH --partition=multix
#SBATCH --nodelist=phoenix3

export OMP_NUM_THREADS=8
hostname

/users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/destab/1nei/activeIsing

