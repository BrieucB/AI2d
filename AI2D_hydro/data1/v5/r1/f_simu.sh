#!/bin/bash
#SBATCH --job-name=PDEv5
#SBATCH -t 13-00:00:00
#SBATCH -n 1
#SBATCH --partition=normal
hostname

srun nice -n 19 python3 /users/invites/benvegnen/Thesis/AI2D/AI2D_hydro/finite_diff_AI2d.py 2000 400 1.5 5 1 2 1 0.001 500 50 2

