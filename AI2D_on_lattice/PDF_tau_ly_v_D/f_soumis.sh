#!/bin/bash

make clean ; make

mkdir data_rev_v0
cd data_rev_v0

for D in 0.5 0.7 1.0
do
    mkdir D$D
    cd D$D

    for ly in 1 2 3 4 5 6 7 8
    do

	mkdir ly${ly}
	cd ly${ly}

	for i in `seq 0 2`
	do
	    mkdir i$i
	    cd i$i
	    
	    cat <<EOF > f_input.dat
tgap = 15000000 tmax = 1000000000 rho0 = 3 lx = 100 ly = ${ly} w0 = 1 beta = 2 v = 1 D = $D

EOF

	    cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=D$D_ly${ly}
#SBATCH -t 14-00:00:00
#SBATCH -n 1
#SBATCH --partition=normal
hostname

/users/invites/benvegnen/AI2D/AI2D_on_lattice/PDF_tau_ly_v_D/activeIsing

EOF
	    chmod u+x f_simu.sh
	    sbatch f_simu.sh
	    cd ..

	done
	cd ..

    done
    cd ..

done
