#!/bin/bash
#SBATCH --job-name=des0.20
#SBATCH -t 7-00:00:00
#SBATCH -n 8
#SBATCH --partition=multi
#SBATCH --nodelist=phoenix0
export OMP_NUM_THREADS=8
hostname

/users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/critical_D/9nei/activeIsing

