#!/bin/bash
#SBATCH --job-name=f_D0.31
#SBATCH -t 7-00:00:00
#SBATCH -n 16
#SBATCH --partition=multi
#SBATCH --nodelist=phoenix3,phoenix4
hostname
export OMP_NUM_THREADS=16

/users/invites/benvegnen/Thesis/AI2D/AI2D_off_lattice/diag_D_rho0/activeIsing

