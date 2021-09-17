#!/bin/bash

make clean ; make

mkdir MSD_v1
cd MSD_v1

for v in 0.1 0.5 1.0 2.0

do
    mkdir v$v
    cd v$v

    cat <<EOF > f_input.dat
tgap = 10000 tmax = 5000 rho0 = 3 lx = 400 ly = 100 w0 = 1 beta = 2 v = $v D = 1
EOF

    cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=b2D1v${v}
#SBATCH -t 7-00:00:00
#SBATCH -n 8
#SBATCH --partition=multi
#SBATCH --nodelist=phoenix
export OMP_NUM_THREADS=8
hostname

/users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/MSD/pupdate/activeIsing

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done
