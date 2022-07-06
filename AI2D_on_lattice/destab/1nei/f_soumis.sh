#!/bin/bash

make clean ; make

mkdir destab_ly
cd destab_ly

for ly in 1 2 3 4 5 6 7

do
    mkdir ly$ly
    cd ly$ly

    cat <<EOF > f_input.dat
tgap = 1000000000 tmax = 2000000 rho0 = 3 lx = 400 ly = ${ly} w0 = 1 beta = 2 v = 1 D = 0.5 phi = 50 rhol = 7 rhog = 1 
EOF

    cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=1n_ly
#SBATCH -t 7-00:00:00
#SBATCH -n 8
#SBATCH --partition=multix
#SBATCH --nodelist=phoenix3

export OMP_NUM_THREADS=8
hostname

/users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/destab/1nei/activeIsing

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done



