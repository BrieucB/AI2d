#!/usr/bin/bash
#SBATCH --job-name=pde10.8
#SBATCH -t 6-00:00:00
#SBATCH --partition=multix96
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=32
hostname

mydir0=$(pwd); # save the current dir name in a variable 

echo ${mydir0};

mydir=${mydir0:7};

mkdir -p /home/$mydir # create a working directory specific for the current job on the host machine

cp $mydir0/f_input.dat /home/$mydir/  # copy the input file in the working dir

cp /users/invites/benvegnen/Thesis/AI2D/AI2D_hydro/PDEs_2d/*.c /home/$mydir/  # copy the *.c in the working dir

cp /users/invites/benvegnen/Thesis/AI2D/AI2D_hydro/PDEs_2d/*.h /home/$mydir/  # copy the *.h in the working dir

cp /users/invites/benvegnen/Thesis/AI2D/AI2D_hydro/PDEs_2d/Makefile /home/$mydir/  # copy the Makefile in the working dir

cd /home/$mydir # go to the working dir

make clean ; make

export OMP_NUM_THREADS=32

srun nice -n 19 ./solve_PDEs_2d

cp *.dat  $mydir0 # copy the result files in your own dir
