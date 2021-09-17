#!/bin/bash
#SBATCH --job-name=td9_0.28
#SBATCH -t 7-00:00:00
#SBATCH -n 8
#SBATCH --partition=multi


export OMP_NUM_THREADS=8
hostname

/users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/critical_D/9nei/activeIsing

