#!/bin/bash
#SBATCH --job-name=PDEv4
#SBATCH -t 13-00:00:00
#SBATCH -n 1
#SBATCH --partition=normal
hostname

srun nice -n 19 python3 /users/invites/benvegnen/Thesis/AI2D/AI2D_hydro/finite_diff_AI2d.py 2000 400 1.5 4 1 2 0 0.001 500 50 2

