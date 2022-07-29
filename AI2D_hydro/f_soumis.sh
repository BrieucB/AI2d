#!/bin/bash

make clean ; make

mkdir data1
cd data1

for v in 1 2 3 4 5
do
    mkdir v$v
    cd v$v
    
    for r in 0 1
    do
	mkdir r$r
	cd r$r

	cat <<EOF > f_simu.sh
#!/bin/bash
#SBATCH --job-name=PDEv$v
#SBATCH -t 13-00:00:00
#SBATCH -n 1
#SBATCH --partition=normal
hostname

srun nice -n 19 python3 /users/invites/benvegnen/Thesis/AI2D/AI2D_hydro/finite_diff_AI2d.py 2000 400 1.5 $v 1 2 $r 0.001 500 50 2

EOF
	chmod u+x f_simu.sh
	sbatch f_simu.sh
	cd ..
    done
    cd ..
done


