#!/bin/bash

make clean ; make

mkdir runs
cd runs

for rho0 in 0.5 1 2 5
do 
    mkdir rho0${rho0}
    cd rho0${rho0}

    mkdir multi
    cd multi
    
    cat <<EOF > f_input.dat
tgap = 10000 tmax = 10000 rho0 = ${rho0} lx = 400 ly = 100 r0 = 1 v0 = 1 D = 1 beta = 2 w0 = 1
	
EOF

    cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=bm_m_${rho0}
#SBATCH -t 7-00:00:00
#SBATCH -n 16
#SBATCH --partition=multi
hostname
export OMP_NUM_THREADS=16

time /users/invites/benvegnen/AI2D/AI2D_off_lattice/benchmark_omp/activeIsing >> time.dat

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh

    cd ..
    mkdir normal
    cd normal

    cat <<EOF > f_input.dat
tgap = 10000 tmax = 10000 rho0 = ${rho0} lx = 400 ly = 100 r0 = 1 v0 = 1 D = 1 beta = 2 w0 = 1
	
EOF
    
    cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=bm_n_${rho0}
#SBATCH -t 7-00:00:00
#SBATCH -n 1
#SBATCH --partition=normal
hostname

time /users/invites/benvegnen/AI2D/AI2D_off_lattice/benchmark_omp/activeIsing >> time.dat

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh

    cd ..
    cd ..
done
