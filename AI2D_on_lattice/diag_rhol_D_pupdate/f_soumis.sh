#!/bin/bash

make clean ; make

mkdir avg_prof_D_v2
cd avg_prof_D_v2

for D in 0.1 0.2 0.3 0.4 0.6 0.8 1.0 1.2 1.4

do
    mkdir D$D
    cd D$D

    cat <<EOF > f_input.dat
tgap = 1000000 tmax = 15000 rho0 = 3 lx = 400 ly = 100 w0 = 1 beta = 2 v = 1 D = $D 
EOF

    cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=bdsD${D}
#SBATCH -t 7-00:00:00
#SBATCH -n 8
#SBATCH --partition=multi
#SBATCH --nodelist=phoenix2
export OMP_NUM_THREADS=8
hostname

/users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/diag_rhol_D_pupdate/activeIsing

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done
