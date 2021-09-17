#!/bin/bash
#SBATCH --job-name=st_D0.4
#SBATCH -t 14-00:00:00
#SBATCH -n 16
#SBATCH --partition=multi
#SBATCH --nodelist=phoenix3,phoenix4
hostname
export OMP_NUM_THREADS=16

/users/invites/benvegnen/AI2D/AI2D_off_lattice/destab_bands/activeIsing

