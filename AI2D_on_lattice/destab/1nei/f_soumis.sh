
#!/bin/bash

make clean ; make

mkdir destab_D
cd destab_D

for D in 0.06 0.08 0.1 0.12 0.14 0.16 0.18

do
    mkdir D$D
    cd D$D

    cat <<EOF > f_input.dat
tgap = 1000000000 tmax = 2000000 rho0 = 3 lx = 200 ly = 100 w0 = 1 beta = 2 D = $D eps = 0.8 phi = 50 rhol = 7 rhog = 1 
EOF

    cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=td9_$D
#SBATCH -t 7-00:00:00
#SBATCH -n 8
#SBATCH --partition=multix

export OMP_NUM_THREADS=8
hostname

/users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/critical_D/9nei/activeIsing

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done


#SBATCH --nodelist=phoenix0
