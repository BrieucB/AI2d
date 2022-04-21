#!/bin/bash
#SBATCH --job-name=des2.140
#SBATCH -t 7-00:00:00
#SBATCH -n 8
#SBATCH --partition=multix

export OMP_NUM_THREADS=8
hostname

srun nice -n 19 /users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/destab_liqp_blob/activeIsing

