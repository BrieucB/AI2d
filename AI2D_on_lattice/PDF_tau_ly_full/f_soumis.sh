#!/bin/bash

make clean ; make

mkdir film_ly
cd film_ly

for ly in 1 2 3 4 5 6 7 

do
    mkdir ly${ly}
    cd ly${ly}

    cat <<EOF > f_input.dat
tgap = 5 tmax = 1000000 rho0 = 3 lx = 100 ly = ${ly} w0 = 1 beta = 2 v = 1 D = 1
EOF

    cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=D.5l${ly}
#SBATCH -t 7-00:00:00
#SBATCH -n 8
#SBATCH --partition=multix
#SBATCH --nodelist=phoenix2

export OMP_NUM_THREADS=8
hostname

/users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/PDF_tau_ly_full/activeIsing

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done


