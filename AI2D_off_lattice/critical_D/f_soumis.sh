#!/bin/bash

make clean ; make

mkdir runs_v3
cd runs_v3

for D in 0.005 0.01 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18 0.2 0.22 0.28 0.3 0.32 0.34 0.36 0.38 0.4 
do
    mkdir D$D
    cd D$D

    cat <<EOF > f_input.dat
tgap = 1 tmax = 10000 rho0 = 3 lx = 200 ly = 100 r0 = 1 v0 = 1 D = $D beta = 2 w0 = 1
	
EOF

    cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=tdoff_${D}
#SBATCH -t 13-00:00:00
#SBATCH -n 16
#SBATCH --partition=normalx
hostname
export OMP_NUM_THREADS=16

srun nice -n 19 /users/invites/benvegnen/Thesis/AI2D/AI2D_off_lattice/critical_D/activeIsing

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done
