#!/usr/bin/bash

make clean ; make

mkdir destab_film_D0.1
cd destab_film_D0.1

cp ../f_soumis.sh .

for i in {0..7} 

do
    mkdir i$i
    cd i$i

    cat <<EOF > f_input.dat
tgap = 1 tmax = 50000000 rho0 = 8 lx = 100 ly = 100 w0 = 1 beta = 2 v = 1 D = 0.1
EOF

    cat <<EOF > f_simu.sh
#!/usr/bin/bash
#SBATCH --job-name=fdD.1
#SBATCH -t 7-00:00:00
#SBATCH --partition=multix64
#SBATCH --nodelist=phoenix8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=8
hostname

mydir0=\$(pwd); # save the current dir name in a variable 

echo \${mydir0};

mydir=\${mydir0:7}

mkdir -p /home/\$mydir # create a working directory specific for the current job on the host machine

cp /users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/spontaneous_destab_liqphase/1box_film_destab/*.c /home/\$mydir/  # copy the input file in the working dir

cp \$mydir0/f_input.dat /home/\$mydir/

cp /users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/spontaneous_destab_liqphase/1box_film_destab/*.h /home/\$mydir/  # copy the executable in the working dir

cp /users/invites/benvegnen/Thesis/AI2D/AI2D_on_lattice/spontaneous_destab_liqphase/1box_film_destab/Makefile /home/\$mydir/

cd /home/\$mydir # go in the in the working dir

make clean ; make

export OMP_NUM_THREADS=8
time srun nice -n 19 ./activeIsing

cp f_profiles*  \$mydir0 # copy the result files in your own dir
cp f_td.dat  \$mydir0 # copy the result files in your own dir
cp f_mag.dat  \$mydir0 # copy the result files in your own dir

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done

