#!/usr/bin/bash

make clean ; make

#rm -rf asympt_prof_varying_beta

mkdir asympt_prof_varying_beta_v4
cd asympt_prof_varying_beta_v4

cp ../f_soumis.sh .

for beta in 1.1 #1.4 1.8 #1.2 1.1 #2.5 

do
    mkdir beta$beta
    cd beta$beta

    cat <<EOF > f_input.dat
tgap = 5000 tmax = 40010 dt = 0.05 lx = 22000 ly = 11000 ds = 1 rhol = 1 beta = $beta v = 1 D = 1 gamma = 1 rhof = 10
EOF

    cat <<EOF > f_simu.sh
#!/usr/bin/bash
#SBATCH --job-name=pde${beta}par
#SBATCH -t 7-00:00:00
#SBATCH --partition=multix96
#SBATCH --nodelist=phoenix6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=64
hostname

mydir0=\$(pwd); # save the current dir name in a variable 

echo \${mydir0};

mydir=\${mydir0:7};

mkdir -p /home/\$mydir # create a working directory specific for the current job on the host machine

cp \$mydir0/f_input.dat /home/\$mydir/  # copy the input file in the working dir

cp /users/invites/benvegnen/Thesis/AI2D/AI2D_hydro/PDEs_2d/solve_PDEs_2d /home/\$mydir/  # copy the executable in the working dir

cd /home/\$mydir # go to the working dir

export OMP_NUM_THREADS=64

srun nice -n 19 ./solve_PDEs_2d

cp *.dat  \$mydir0 # copy the result files in your own dir
EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done
