#!/bin/bash

make clean ; make

mkdir data_varying_r

cp f_soumis.sh data_varying_r

cd data_varying_r

for r in 0.0 0.2 0.4 0.6 0.8 1.0

do
    mkdir r$r
    cd r$r
    cat <<EOF > f_input.dat
tgap = 20 tmax = 2000 dt = 0.01 lx = 10000 dx = 0.5 rhol = 1 beta = 2 v = 1 D = 1 r = $r rhoc = 3
EOF

    cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=pde1d
#SBATCH -t 13-00:00:00
#SBATCH -n 1
#SBATCH --nice=19
#SBATCH --partition=normal

hostname

srun nice -n 19 /users/invites/benvegnen/Thesis/AI1D/PDEs_1d/solve_PDEs_1d

EOF
    chmod u+x f_simu.sh
    sbatch f_simu.sh
    cd ..
done
