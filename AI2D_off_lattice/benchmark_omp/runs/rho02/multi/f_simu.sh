#!/bin/bash
#SBATCH --job-name=bm_m_2
#SBATCH -t 7-00:00:00
#SBATCH -n 16
#SBATCH --partition=multi
hostname
export OMP_NUM_THREADS=16

time /users/invites/benvegnen/AI2D/AI2D_off_lattice/benchmark_omp/activeIsing >> time.dat

