#!/bin/bash

make clean ; make

mkdir avg_prof_D_v2
cd avg_prof_D_v2

for D in 0.1 0.2 0.3 0.4 0.6 0.8 1.0 1.2 1.4 # 1.6 1.8 2.0 2.2 2.4 2.6

do
    mkdir D$D
    cd D$D

    cat <<EOF > f_input.dat
tgap = 10 tmax = 15000 rho0 = 3 lx = 400 ly = 100 w0 = 1 beta = 2 v = 1 D = $D 
EOF

    cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=sb2v1D${D}
#SBATCH -t 7-00:00:00
#SBATCH -n 1
#SBATCH --partition=normal
hostname

/users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/diag_rhol_D_supdate/activeIsing

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done
