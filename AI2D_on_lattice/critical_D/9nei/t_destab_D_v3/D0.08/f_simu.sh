#!/bin/bash
#SBATCH --job-name=td9_0.08
#SBATCH -t 13-00:00:00
#SBATCH -n 8
#SBATCH --partition=normalx


export OMP_NUM_THREADS=8
hostname

/users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/critical_D/9nei/activeIsing

