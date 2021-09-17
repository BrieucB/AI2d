#!/bin/bash

make clean ; make

mkdir avg_prof_D_v3
cd avg_prof_D_v3

for D in 0.31 0.63 0.94 1.26 1.88 2.51 3.14 3.77 #0.1 0.2 0.3 0.4 0.6 0.8 1.0 1.2

do
    mkdir D$D
    cd D$D

    cat <<EOF > f_input.dat
tgap = 100 tmax = 5000 rho0 = 0.95 lx = 708 ly = 177 r0 = 1 v0 = 1.77 D = $D beta = 2 w0 = 1
	
EOF

    cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=f_D${D}
#SBATCH -t 7-00:00:00
#SBATCH -n 16
#SBATCH --partition=multi
#SBATCH --nodelist=phoenix3,phoenix4
hostname
export OMP_NUM_THREADS=16

/users/invites/benvegnen/Thesis/AI2D/AI2D_off_lattice/diag_D_rho0/activeIsing

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done
