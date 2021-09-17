#!/bin/bash
#SBATCH --job-name=bm_n_2
#SBATCH -t 7-00:00:00
#SBATCH -n 1
#SBATCH --partition=normal
hostname

time /users/invites/benvegnen/AI2D/AI2D_off_lattice/benchmark_omp/activeIsing >> time.dat

