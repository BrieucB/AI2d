#!/bin/bash

make clean ; make

mkdir instances_beta
cd instances_beta

for i in "1.05 29" #"1.1 21.53585" "1.2 11.361049999999999" "1.4 6.60145" "1.6 4.9011949999999995" "1.9 3.522575" "2.3 2.4504455" "2.5 2.1005" "3.0 2"
do
    set -- $i
    
    mkdir beta$1
    cd beta$1
    
    cat <<EOF > f_input.dat
tgap = 1500 tmax = 15000 rho0 = $2 lx = 400 ly = 400 w0 = 1 beta = $1 eps = 0.9 D = 1

EOF

    cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=bands_$1
#SBATCH -t 7-00:00:00
#SBATCH -n 1
#SBATCH --partition=normal
hostname

/users/invites/benvegnen/AI2D/diag_T_rho/activeIsing

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done
