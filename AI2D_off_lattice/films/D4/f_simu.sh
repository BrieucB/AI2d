#!/bin/bash
#SBATCH --job-name=f_D4
#SBATCH -t 7-00:00:00
#SBATCH -n 16
#SBATCH --partition=normal
hostname
export OMP_NUM_THREADS=16

/users/invites/benvegnen/AI2D/AI2D_off_lattice/activeIsing

