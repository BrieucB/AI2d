#!/bin/bash

make clean ; make

mkdir MSD_v2
cd MSD_v2

for v in 0.1 0.5 1.0 2.0

do
    mkdir v$v
    cd v$v

    cat <<EOF > f_input.dat
tgap = 100000 tmax = 5000 rho0 = 3 lx = 400 ly = 100 w0 = 1 beta = 2 v = $v D = 0 
EOF

    cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=f_D${D}
#SBATCH -t 7-00:00:00
#SBATCH -n 1
#SBATCH --partition=normal
hostname

/users/invites/benvegnen/AI2D/AI2D_on_lattice/MSD/supdate/activeIsing

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done
