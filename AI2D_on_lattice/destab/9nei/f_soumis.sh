
#!/bin/bash

make clean ; make

mkdir t_destab_D_v2
cd t_destab_D_v2

for D in 

do
    mkdir D$D
    cd D$D

    cat <<EOF > f_input.dat
tgap = 1 tmax = 2000000 rho0 = 3 lx = 200 ly = 100 w0 = 1 beta = 2 v = 1 D = $D phi = 50 rhol = 7 rhog = 1 
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


#SBATCH --nodelist=phoenix0
