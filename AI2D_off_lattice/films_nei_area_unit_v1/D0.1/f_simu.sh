#!/bin/bash
#SBATCH --job-name=f1_D0.1
#SBATCH -t 7-00:00:00
#SBATCH -n 16
#SBATCH --partition=multi
#SBATCH --nodelist=phoenix3
hostname
export OMP_NUM_THREADS=16

/users/invites/benvegnen/Thesis/AI2D/AI2D_off_lattice/activeIsing

