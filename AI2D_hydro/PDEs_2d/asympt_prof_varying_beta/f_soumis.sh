#!/usr/bin/bash

make clean ; make

rm -rf asympt_prof_varying_beta

mkdir asympt_prof_varying_beta
cd asympt_prof_varying_beta

cp ../f_soumis.sh .

for beta in 4.0 2.5 1.8 1.4 1.2 1.01

do
    mkdir beta$beta
    cd beta$beta

    cat <<EOF > f_input.dat
tgap = 1000 tmax = 10000 dt = 0.05 lx = 3000 ly = 2000 ds = 1 rhol = 1 beta = $beta v = 1 D = 1 gamma = 1 rhof = 10
EOF

    cat <<EOF > f_simu.sh
#!/usr/bin/bash
#SBATCH --job-name=pde$beta
#SBATCH -t 7-00:00:00
#SBATCH -n 1
#SBATCH --partition=normal

hostname

mydir0=\$(pwd); # save the current dir name in a variable 

echo \${mydir0};

mydir=\${mydir0:7};

mkdir -p /home/\$mydir # create a working directory specific for the current job on the host machine

cp \$mydir0/f_input.dat /home/\$mydir/  # copy the input file in the working dir

cp /users/invites/benvegnen/Thesis/AI2D/AI2D_hydro/PDEs_2d/solve_PDEs_2d /home/\$mydir/  # copy the executable in the working dir

cd /home/\$mydir # go to the working dir

srun nice -n 19 ./solve_PDEs_2d

cp *.dat  \$mydir0 # copy the result files in your own dir
EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done
