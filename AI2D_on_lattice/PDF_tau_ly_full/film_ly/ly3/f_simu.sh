#!/bin/bash
#SBATCH --job-name=D.5l3
#SBATCH -t 7-00:00:00
#SBATCH -n 8
#SBATCH --partition=multix
#SBATCH --nodelist=phoenix2

export OMP_NUM_THREADS=8
hostname

/users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/PDF_tau_ly_full/activeIsing

