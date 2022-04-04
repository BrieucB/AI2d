
#!/bin/bash

make clean ; make

mkdir t_destab_D_avg_m
cd t_destab_D_avg_m

for D in 0.26 #0.005 0.01 0.22 0.24 0.26 0.28 0.3 0.32 0.34 0.36 # 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18 0.2 #0.32 0.34 0.36 0.38 0.4 #0.24 0.26 0.28 0.3 #0.01 0.005 0.22 # # 

do
    mkdir D$D
    cd D$D

    cat <<EOF > f_input.dat
tgap = 1000000000 tmax = 200000000 rho0 = 3 lx = 200 ly = 100 w0 = 1 beta = 2 v = 1 D = $D phi = 50 rhol = 7 rhog = 1 
EOF

    cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=td9_$D
#SBATCH -t 7-00:00:00
#SBATCH -n 8
#SBATCH --partition=multi


export OMP_NUM_THREADS=8
hostname

/users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/critical_D/9nei/activeIsing

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done




#SBATCH --nodelist=phoenix1
