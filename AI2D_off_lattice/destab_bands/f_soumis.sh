#!/bin/bash

make clean ; make

mkdir runs_v1
cd runs_v1

for D in 0.05 0.08 0.1 0.2 0.3 0.4 0.5 0.6 0.7
do
    mkdir D$D
    cd D$D

    cat <<EOF > f_input.dat
tgap = 1 tmax = 10000 rho0 = 1 lx = 400 ly = 100 r0 = 0.564 v0 = 1 D = $D beta = 2 w0 = 1
	
EOF

    cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=st_D${D}
#SBATCH -t 14-00:00:00
#SBATCH -n 16
#SBATCH --partition=multi
#SBATCH --nodelist=phoenix3,phoenix4
hostname
export OMP_NUM_THREADS=16

/users/invites/benvegnen/AI2D/AI2D_off_lattice/destab_bands/activeIsing

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done
