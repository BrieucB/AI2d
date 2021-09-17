#!/bin/bash

make clean ; make

mkdir films_nei_area_unit_v1
cd films_nei_area_unit_v1

for D in 0.4 0.6 0.8 1.0 1.5
do
    mkdir D$D
    cd D$D

    cat <<EOF > f_input.dat
tgap = 10 tmax = 5000 rho0 = 0.95 lx = 708 ly = 177 r0 = 1 v0 = 1.77 D = $D beta = 2 w0 = 1
	
EOF

    cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=f1_D${D}
#SBATCH -t 7-00:00:00
#SBATCH -n 16
#SBATCH --partition=multi
#SBATCH --nodelist=phoenix
hostname
export OMP_NUM_THREADS=16

/users/invites/benvegnen/Thesis/AI2D/AI2D_off_lattice/activeIsing

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done
