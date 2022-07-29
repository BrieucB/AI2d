#!/bin/bash
#SBATCH --job-name=PDEv1
#SBATCH -t 13-00:00:00
#SBATCH -n 1
#SBATCH --partition=normal
hostname

srun nice -n 19 python3 /users/invites/benvegnen/Thesis/AI2D/AI2D_hydro/finite_diff_AI2d.py 10000 2000 1.5 1 1 2 1 0.01 2200 300 2

