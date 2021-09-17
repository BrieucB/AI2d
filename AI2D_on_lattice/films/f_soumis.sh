#!/bin/bash

make clean ; make

mkdir instances_Dy
cd instances_Dy

for Dy in 0.025 0.030 0.035 0.040 0.045 
do
    mkdir Dy${Dy}
    cd Dy${Dy}
    
    cat <<EOF > f_input.dat
tgap = 10 tmax = 10000 rho0 = 3 lx = 400 ly = 50 w0 = 1 beta = 2 eps = 0.9 Dx = 1 Dy = ${Dy}

EOF

    cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=f_D${Dy}
#SBATCH -t 7-00:00:00
#SBATCH -n 1
#SBATCH --partition=normal
hostname

/users/invites/benvegnen/AI2D/films/activeIsing

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done
