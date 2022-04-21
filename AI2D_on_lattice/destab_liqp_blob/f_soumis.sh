#!/bin/bash

make clean ; make

mkdir data_beta_h 
cd data_beta_h
 
for beta in 2.7 2.9 3.1 3.3 3.5 # 1.2 1.5 1.7 1.9 2.1 2.3 2.5 
do
    mkdir beta$beta
    cd beta$beta

    for h in {10..150..10}
    do

	mkdir h$h
	cd h$h
	
	cat <<EOF > f_input.dat
tgap = 100 tmax = 100 rho0 = 5 lx = 100 ly = 20 w0 = 1 beta = $beta v = 1 D = 1 h0 = $h
EOF

	cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=des${beta}$h
#SBATCH -t 7-00:00:00
#SBATCH -n 8
#SBATCH --partition=multix

export OMP_NUM_THREADS=8
hostname

srun nice -n 19 /users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/destab_liqp_blob/activeIsing

EOF
	chmod u+x f_simu.sh
	sbatch f_simu.sh
	cd ..
    done
    cd ..
done
