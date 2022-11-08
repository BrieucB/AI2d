#!/usr/bin/bash

make clean ; make

mkdir critical_rhof
cd critical_rhof

ncpu=12

cp ../f_soumis.sh .

for H in 1.0 1.2 1.5 1.7 2.0 2.2 2.5 2.7 3.0 3.2 3.5 3.7 4.0

do
    mkdir H$H
    cd H$H

    cat <<EOF > f_input.dat
ncpu = $ncpu nby_box = 4 tgap = 1 tmax = 41 dt = 0.05 lx = 90 ly = 90 ds = 0.5 rhol = 1 beta = 2 v = 1 D = 1 gamma = 1 rhof = $H
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
