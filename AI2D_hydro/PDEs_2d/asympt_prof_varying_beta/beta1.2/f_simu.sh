#!/usr/bin/bash
#SBATCH --job-name=pde1.2
#SBATCH -t 7-00:00:00
#SBATCH -n 1
#SBATCH --partition=normal

hostname

mydir0=$(pwd); # save the current dir name in a variable 

echo ${mydir0};

mydir=${mydir0:7};

mkdir -p /home/$mydir # create a working directory specific for the current job on the host machine

cp $mydir0/f_input.dat /home/$mydir/  # copy the input file in the working dir

cp /users/invites/benvegnen/Thesis/AI2D/AI2D_hydro/PDEs_2d/solve_PDEs_2d /home/$mydir/  # copy the executable in the working dir

cd /home/$mydir # go to the working dir

srun nice -n 19 ./solve_PDEs_2d

cp *.dat  $mydir0 # copy the result files in your own dir
