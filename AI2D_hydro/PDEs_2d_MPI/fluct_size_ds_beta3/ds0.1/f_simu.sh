#!/usr/bin/bash
#SBATCH --job-name=pde2dx0.1
#SBATCH -t 7-00:00:00
#SBATCH --partition=multix96,multix64,multix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=24
hostname

mydir0=$(pwd); # save the current dir name in a variable 

echo ${mydir0};

mydir=${mydir0:7};

mkdir -p /home/$mydir # create a working directory specific for the current job on the host machine

cp $mydir0/f_input.dat /home/$mydir/  # copy the input file in the working dir

cp /users/invites/benvegnen/Thesis/AI2D/AI2D_hydro/PDEs_2d_MPI/*.c /home/$mydir/  # copy the .c files in the working dir

cp /users/invites/benvegnen/Thesis/AI2D/AI2D_hydro/PDEs_2d_MPI/*.h /home/$mydir/  # copy the .h files in the working dir

cp /users/invites/benvegnen/Thesis/AI2D/AI2D_hydro/PDEs_2d_MPI/Makefile /home/$mydir/  # copy the Makefile in the working dir
 
cd /home/$mydir # go to the working dir

make clean ; make

mpiexec --oversubscribe -np  24 ./solve_PDEs_2d #srun nice -n 19 

cp *.dat  $mydir0 # copy the result files in your own dir
