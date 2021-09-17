#!/bin/bash

make clean ; make

mkdir MSD_v2
cd MSD_v2

for v in 0.1 0.5 1.0 2.0

do
    mkdir v$v
    cd v$v

    cat <<EOF > f_input.dat
tgap = 10000 tmax = 5000 rho0 = 1 lx = 400 ly = 100 r0 = 1 v0 = $v D = 0 beta = 2 w0 = 1
EOF

    cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=MSD_v${v}
#SBATCH -t 7-00:00:00
#SBATCH -n 16
#SBATCH --partition=multi
#SBATCH --nodelist=phoenix3,phoenix2
export OMP_NUM_THREADS=16
hostname

/users/invites/benvegnen/AI2D/AI2D_off_lattice/MSD/activeIsing

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done
