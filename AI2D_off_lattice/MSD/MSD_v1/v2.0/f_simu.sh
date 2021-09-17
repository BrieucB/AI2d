#!/bin/bash
#SBATCH --job-name=MSD_v2.0
#SBATCH -t 7-00:00:00
#SBATCH -n 16
#SBATCH --partition=multi
#SBATCH --nodelist=phoenix3,phoenix2
export OMP_NUM_THREADS=16
hostname

/users/invites/benvegnen/AI2D/AI2D_off_lattice/MSD/activeIsing

