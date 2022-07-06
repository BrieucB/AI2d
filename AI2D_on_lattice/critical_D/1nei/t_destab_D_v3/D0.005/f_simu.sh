#!/bin/bash
#SBATCH --job-name=td1_0.005
#SBATCH -t 13-00:00:00
#SBATCH -n 8
#SBATCH --partition=normalx

export OMP_NUM_THREADS=8
hostname

/users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/critical_D/1nei/activeIsing

