#!/bin/bash

make clean ; make

mkdir avg_prof_D
cd avg_prof_D

for D in 0.1 0.12 0.14 0.16 0.18

do
    mkdir D$D
    cd D$D

    cat <<EOF > f_input.dat
tgap = 1000000 tmax = 15000 rho0 = 5 lx = 400 ly = 100 w0 = 1 beta = 2 eps = 0.8 D = $D 
EOF

    cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=bdsD${D}
#SBATCH -t 7-00:00:00
#SBATCH --partition=multix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=8
hostname

export OMP_NUM_THREADS=8

srun nice -n 19 /users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/diag_rho_D_AIMrates/activeIsing

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done
