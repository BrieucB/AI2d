#!/usr/bin/bash

mkdir droplet_size_ds_2
cd droplet_size_ds_2

ncpu=32

cp ../f_soumis.sh .

for ds in 1.0 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 

do
    mkdir ds$ds
    cd ds$ds

    cat <<EOF > f_input.dat
tgap = 2000 tmax = 2100 dt = 0.01 lx = 1500 ly = 600 ds = $ds rhol = 1 beta = 2.5 v = 1 D = 0.5 gamma = 1 rhof = 10
EOF

    cat <<EOF > f_simu.sh
#!/usr/bin/bash
#SBATCH --job-name=pde1${ds}
#SBATCH -t 6-00:00:00
#SBATCH --partition=multix96
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=$ncpu
hostname

mydir0=\$(pwd); # save the current dir name in a variable 

echo \${mydir0};

mydir=\${mydir0:7};

mkdir -p /home/\$mydir # create a working directory specific for the current job on the host machine

cp \$mydir0/f_input.dat /home/\$mydir/  # copy the input file in the working dir

cp /users/invites/benvegnen/Thesis/AI2D/AI2D_hydro/PDEs_2d/*.c /home/\$mydir/  # copy the *.c in the working dir

cp /users/invites/benvegnen/Thesis/AI2D/AI2D_hydro/PDEs_2d/*.h /home/\$mydir/  # copy the *.h in the working dir

cp /users/invites/benvegnen/Thesis/AI2D/AI2D_hydro/PDEs_2d/Makefile /home/\$mydir/  # copy the Makefile in the working dir

cd /home/\$mydir # go to the working dir

make clean ; make

export OMP_NUM_THREADS=$ncpu

srun nice -n 19 ./solve_PDEs_2d

cp *.dat  \$mydir0 # copy the result files in your own dir
EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done
