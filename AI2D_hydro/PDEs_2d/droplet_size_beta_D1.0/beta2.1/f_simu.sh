#!/usr/bin/bash
#SBATCH --job-name=pde12.1
#SBATCH -t 7-00:00:00
#SBATCH --partition=multix64
#SBATCH --nodelist=phoenix11
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=12
hostname

mydir0=$(pwd); # save the current dir name in a variable 

echo ${mydir0};

mydir=${mydir0:7};

mkdir -p /home/$mydir # create a working directory specific for the current job on the host machine

cp $mydir0/f_input.dat /home/$mydir/  # copy the input file in the working dir

cp /users/invites/benvegnen/Thesis/AI2D/AI2D_hydro/PDEs_2d/solve_PDEs_2d /home/$mydir/  # copy the executable in the working dir

cd /home/$mydir # go to the working dir

export OMP_NUM_THREADS=12

srun nice -n 19 ./solve_PDEs_2d

cp *.dat  $mydir0 # copy the result files in your own dir
