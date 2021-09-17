#!/bin/bash

make clean ; make

mkdir instances_lx100_D0.2
cd instances_lx100_D0.2

for ly in 20 50
do
    mkdir ly${ly}
    cd ly${ly}
    for i in `seq 0 15`
    do
	mkdir i$i
	cd i$i

	cat <<EOF > f_input.dat
tgap = 1000000 tmax = 10000000 rho0 = 2 lx = 100 ly = ${ly} w0 = 1 beta = 2 v = 1 D = 0.2
EOF

	cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=pdf_ly${ly}i$i
#SBATCH -t 7-00:00:00
#SBATCH -n 1
#SBATCH --partition=normal
hostname

/users/invites/benvegnen/AI2D/PDF_tau_ly/activeIsing

EOF
	chmod u+x f_simu.sh
	sbatch f_simu.sh
	cd ..
    done
    cd ..
done
