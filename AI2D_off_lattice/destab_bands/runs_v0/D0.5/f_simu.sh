#!/bin/bash
#SBATCH --job-name=st_D0.5
#SBATCH -t 7-00:00:00
#SBATCH -n 16
#SBATCH --partition=multi
hostname
export OMP_NUM_THREADS=16

/users/invites/benvegnen/AI2D/AI2D_off_lattice/destab_bands/activeIsing

