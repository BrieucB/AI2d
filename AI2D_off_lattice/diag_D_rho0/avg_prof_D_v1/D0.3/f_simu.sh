#!/bin/bash
#SBATCH --job-name=f_D0.3
#SBATCH -t 7-00:00:00
#SBATCH -n 16
#SBATCH --partition=multi
#SBATCH --nodelist=phoenix3,phoenix2
hostname
export OMP_NUM_THREADS=16

/users/invites/benvegnen/AI2D/AI2D_off_lattice/diag_D_rho0/activeIsing

