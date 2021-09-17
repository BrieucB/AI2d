#!/bin/bash

make clean ; make

mkdir data_pics_rho5
cd data_pics_rho5

for D in 0.005 0.008 0.01 0.02 0.05 0.08 #0.1 0.2 0.3 0.4 0.5
 do
    mkdir D$D
    cd D$D

    for ly in 100
    do

	mkdir ly${ly}
	cd ly${ly}
	
	cat <<EOF > f_input.dat
tgap = 10 tmax = 2000 rho0 = 5 lx = 400 ly = ${ly} w0 = 1 beta = 2 v = 1 D = $D

EOF

	cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=D${D}_ly${ly}
#SBATCH -t 14-00:00:00
#SBATCH -n 1
#SBATCH --partition=normal
hostname

/users/invites/benvegnen/AI2D/AI2D_on_lattice/films_v_D/activeIsing

EOF
	chmod u+x f_simu.sh
	sbatch f_simu.sh
	cd ..

    done
    cd ..

done
