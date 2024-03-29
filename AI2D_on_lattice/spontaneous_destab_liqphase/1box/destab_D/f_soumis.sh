#!/usr/bin/bash

make clean ; make

###################################
rm -rf destab_D
###################################

mkdir destab_D
cd destab_D

cp ../f_soumis.sh .

for D in 0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 

do
    mkdir D$D
    cd D$D

    cat <<EOF > f_input.dat
tgap = 1 tmax = 2000000 rho0 = 10 lx = 2000 ly = 50 w0 = 1 beta = 2 v = 1 D = $D
EOF

    cat <<EOF > f_simu.sh
#!/usr/bin/bash
#SBATCH --job-name=td1n
#SBATCH -t 7-00:00:00
#SBATCH -n 8
#SBATCH --partition=normalx

hostname

mydir0=\$(pwd); # save the current dir name in a variable 

echo \${mydir0};

mydir=\${mydir0:7}

mkdir -p /home/\$mydir # create a working directory specific for the current job on the host machine

cp /users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/spontaneous_destab_liqphase/1box/f_input.dat /home/\$mydir/  # copy the input file in the working dir

cp /users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/spontaneous_destab_liqphase/1box/activeIsing /home/\$mydir/  # copy the executable in the working dir

cd /home/\$mydir # go in the in the working dir

export OMP_NUM_THREADS=8
srun nice -n 19 ./activeIsing

cp f_profiles*  \$mydir0 # copy the result files in your own dir
cp f_td.dat  \$mydir0 # copy the result files in your own dir
cp f_mag.dat  \$mydir0 # copy the result files in your own dir

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done

