#!/bin/bash

make clean ; make

mkdir films_D
cd films_D

for D in 0.02 0.03 0.05 0.06 0.13 0.19 0.25 0.32 0.48

do
    mkdir D$D
    cd D$D

    cat <<EOF > f_input.dat
tgap = 10 tmax = 5000 rho0 = 3 lx = 400 ly = 100 w0 = 1 beta = 2 v = 1 D = $D 
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

/users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/films_pupdate/activeIsing

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done
