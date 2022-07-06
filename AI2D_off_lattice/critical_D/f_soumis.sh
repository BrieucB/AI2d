#!/bin/bash

make clean ; make

mkdir runs_v2
cd runs_v2

for D in 0.5 0.6 0.7 #0.05 0.08 0.1 0.2 0.3 0.4
do
    mkdir D$D
    cd D$D

    cat <<EOF > f_input.dat
tgap = 1 tmax = 10000 rho0 = 3 lx = 200 ly = 100 r0 = 1 v0 = 1 D = $D beta = 2 w0 = 1
	
EOF

    cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=st_D${D}
#SBATCH -t 7-00:00:00
#SBATCH -n 16
#SBATCH --partition=multix
hostname
export OMP_NUM_THREADS=16

srun nice -n 19 /users/invites/benvegnen/Thesis/AI2D/AI2D_off_lattice/destab_bands/activeIsing

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done
