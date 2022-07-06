#!/bin/bash
#SBATCH --job-name=1n_0.1
#SBATCH -t 7-00:00:00
#SBATCH -n 8
#SBATCH --partition=multix

export OMP_NUM_THREADS=8
hostname

/users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/destab/1nei/activeIsing

