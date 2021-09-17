#!/bin/bash
#SBATCH --job-name=f_D
#SBATCH -t 7-00:00:00
#SBATCH -n 16
#SBATCH --partition=multi
#SBATCH --nodelist=phoenix3
export OMP_NUM_THREADS=16
hostname

/users/invites/benvegnen/AI2D/AI2D_off_lattice/MSD/activeIsing

