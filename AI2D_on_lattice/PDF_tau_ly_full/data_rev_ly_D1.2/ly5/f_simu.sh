#!/bin/bash
#SBATCH --job-name=r_.7D5
#SBATCH -t 7-00:00:00
#SBATCH -n 8
#SBATCH --partition=multi

export OMP_NUM_THREADS=8
hostname

/users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/PDF_tau_ly_full/activeIsing

