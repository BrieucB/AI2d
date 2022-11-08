#!/usr/bin/bash

make clean ; make

mkdir fluct_size_ds
cd fluct_size_ds

ncpu=32

cp ../f_soumis.sh .

for ds in 0.1 0.5 1.0

do
    mkdir ds$ds
    cd ds$ds

    cat <<EOF > f_input.dat
ncpu = $ncpu nby_box = 4 tgap = 2000 tmax = 2100 dt = 0.01 lx = 2000 ly = 1000 ds = $ds rhol = 1 beta = 2 v = 1 D = 1 gamma = 1 rhof = 10
EOF

    cat <<EOF > f_simu.sh
#!/usr/bin/bash
#SBATCH --job-name=pde2b${beta}
#SBATCH -t 7-00:00:00
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

cp /users/invites/benvegnen/Thesis/AI2D/AI2D_hydro/PDEs_2d_MPI/solve_PDEs_2d /home/\$mydir/  # copy the executable in the working dir

cd /home/\$mydir # go to the working dir

srun nice -n 19 mpirun -oversubscribe -np $ncpu ./solve_PDEs_2d

cp *.dat  \$mydir0 # copy the result files in your own dir
EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done
