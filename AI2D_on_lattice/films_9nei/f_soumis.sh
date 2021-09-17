#!/bin/bash

make clean ; make

mkdir films_D_adim
cd films_D_adim

for D in 0.18 0.27 0.45 0.54 1.17 1.71 2.24 2.88 4.32 #0.02 0.03 0.05 0.06 0.13 0.19 0.25 0.32 0.48

do
    mkdir D$D
    cd D$D

    cat <<EOF > f_input.dat
tgap = 10 tmax = 5000 rho0 = 1 lx = 1200 ly = 300 w0 = 1 beta = 2 v = 3 D = $D 
EOF

    cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=f_D${D}
#SBATCH -t 7-00:00:00
#SBATCH -n 8
#SBATCH --partition=multi
#SBATCH --nodelist=phoenix2
export OMP_NUM_THREADS=8
hostname

/users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/films_9nei/activeIsing

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done
