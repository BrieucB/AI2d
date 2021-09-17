#!/bin/bash

make clean ; make

mkdir instances_beta
cd instances_beta

for beta in 2.5 2.2 2.0 1.8 1.6 1.5 1.4 1.3 1.2 1.1
do
    mkdir beta${beta}
    cd beta${beta}
    
    cat <<EOF > f_input.dat
tgap = 1500 tmax = 15000 rho0 = 3 lx = 400 ly = 400 w0 = 1 beta = ${beta} v = 1 D = 1

EOF

    cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=v_b${beta}
#SBATCH -t 7-00:00:00
#SBATCH -n 1
#SBATCH --partition=normal
hostname

/users/invites/benvegnen/AI2D/band_velocity/activeIsing

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done
