#!/usr/bin/bash

make clean ; make

#rm -rf asympt_prof_varying_beta

mkdir droplet_size_beta
cd droplet_size_beta

ncpu=64

cp ../f_soumis.sh .

for beta in 1.08 #2.6 2.7 2.8 2.9 3.0 #2.1 2.2 2.3 2.4 2.5 #1.6 1.7 1.8 1.9 2.0 #1.1 1.2 1.3 1.4 1.5

do
    mkdir beta$beta
    cd beta$beta

    cat <<EOF > f_input.dat
tgap = 5000 tmax = 5100 dt = 0.05 lx = 5500 ly = 2200 ds = 1 rhol = 1 beta = $beta v = 1 D = 1 gamma = 1 rhof = 10
EOF

    cat <<EOF > f_simu.sh
#!/usr/bin/bash
#SBATCH --job-name=pde1${beta}
#SBATCH -t 7-00:00:00
#SBATCH --partition=multix64
#SBATCH --nodelist=phoenix11
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

cp /users/invites/benvegnen/Thesis/AI2D/AI2D_hydro/PDEs_2d/solve_PDEs_2d /home/\$mydir/  # copy the executable in the working dir

cd /home/\$mydir # go to the working dir

export OMP_NUM_THREADS=$ncpu

srun nice -n 19 ./solve_PDEs_2d

cp *.dat  \$mydir0 # copy the result files in your own dir
EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done
