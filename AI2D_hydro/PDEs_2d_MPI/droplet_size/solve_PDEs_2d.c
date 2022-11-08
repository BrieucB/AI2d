#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#include "allocation.h"
#include "exit_if.h"
#include "erreurs.h"
#include "tools.h"

int main(int argc, char** argv)
{
  // MPI INIT
  int rank, numtasks;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  // FILES
  FILE *f_input, *f_rho, *f_m;
  f_input=fopen("f_input.dat", "r");

  // PARAMETERS
  int lx, ly;
  int tmax;
  double tgap;
  
  double dt;
  double ds;
  double beta;
  double D;
  double v;
  double rhol;
  double gamma;
  double rhof;

  int ncpu;
  int nby_box;
  
  fscanf(f_input, "ncpu = %d nby_box = %d tgap = %lg tmax = %d dt = %lg lx = %d ly = %d ds = %lg rhol = %lg beta = %lg v = %lg D = %lg gamma = %lg rhof = %lg", &ncpu, &nby_box, &tgap, &tmax, &dt, &lx, &ly, &ds, &rhol, &beta, &v, &D, &gamma, &rhof);

  /* printf("ncpu = %d tgap = %lg tmax = %d dt = %lg lx = %d ly = %d ds = %lg rhol = %lg beta = %lg v = %lg D = %lg gamma = %lg rhof = %lg\n", ncpu, tgap, tmax, dt, lx, ly, ds, rhol, beta, v, D, gamma, rhof); */
  

  if(ncpu!=numtasks)
    {
      printf("Wrong number of cores allocated. numtasks = %d\n", numtasks);
      exit(0);
    }

  // MPI BOXES
  int nbx_box = (int) ncpu/nby_box; // ncpu is a multiple of nby_box
  
  /* SPACE-TIME DISCRETIZATION */
  double Nt= floor(tmax/dt);
  int Nx= (int) floor(lx/ds);
  Nx = Nx - Nx%nbx_box;
  int Ny= (int) floor(ly/ds);
  Ny = Ny - Ny%nby_box;
  int Ny2= (int) floor(Ny/2);
  int Nout=(int) floor(tgap/dt);
  int clap = Nout;
  double c_diff=D*(dt/(ds*ds));
  double c_adv=(v/2.)*(dt/ds);

  // => nby_box*nbx_box = ncpu
  
  int Nx_box = (int) Nx/nbx_box;
  int Ny_box = (int) Ny/nby_box;

  /* Droplet's tail cut */
  int rank_x, x0_shifted, box_x, x0_rel, x01, box_x1, x01_rel;
  int shift=(int) floor(10/ds);
  
  int box0, old_box0;
  int x0, old_x0, y0;

  /* Increments */
  int y;
  
  /* if(rank==0) */
  /*   printf("Nx = %d, Ny = %d, Nx_box = %d, Ny_box = %d\n", Nx, Ny, Nx_box, Ny_box); */
  
  double **rho_loc = dmatrix((long) Nx_box +2, (long) Ny_box +2);
  double **m_loc = dmatrix((long) Nx_box +2, (long) Ny_box +2);
  
  double **rho_loc_old = dmatrix((long) Nx_box +2, (long) Ny_box +2);
  double **m_loc_old = dmatrix((long) Nx_box +2, (long) Ny_box +2);

  /* Matrices containing the whole lattice */ 
  double **rho_tot;
  double **m_tot;

  /* Matrices to broadcast the entire *_tot matrices on the local *_loc cpus */
  double **m_rec;
  double **rho_rec;

  /* Matrices to gather the local *_loc cpus on the global *_tot matrices (cpu 0) */
  double **m_rec0 = dmatrix((long) Nx_box +2, (long) Ny_box +2);
  double **rho_rec0 = dmatrix((long) Nx_box +2, (long) Ny_box +2);;

  /* Vector to send & receive the horizontal borders */
  MPI_Datatype sent_hborder;
  double *rec_hborder =  (double *)malloc((Nx_box+2)*sizeof(double));
  
  int *frame = (int *)malloc(4*nbx_box*nby_box*sizeof(int));
  init_frame(nbx_box, nby_box, frame);
  
  //Determine ml = rhol*tanh(beta*phil)
  double phil = 0.5;
  for(int p=0 ; p<100 ; p++)
    {
      phil=tanh(beta*phil);
    }
  double ml=rhol*phil;

  /* Allocate memory for node 0 */
  if(rank==0)
    {
      rho_tot = dmatrix((long) Nx, (long) Ny);
      m_tot = dmatrix((long) Nx, (long) Ny);
    }

  /* m_rec = dmatrix((long) Nx_box+2, (long) Ny_box+2); */
  /* rho_rec = dmatrix((long) Nx_box+2, (long) Ny_box+2); */

  m_rec = dmatrix((long) Nx_box, (long) Ny_box);
  rho_rec = dmatrix((long) Nx_box, (long) Ny_box);
  
  /* Initialize fields */
  initFields(rank, rho_tot, m_tot, m_rec, Nx, Ny, Ny2, Nx_box, Ny_box, ds, rhof, rhol,ml);

  if(rank==0)
    printf("Init done\n");

  /* Scatter fields among the cpus */
  scatterTotalMatrices(rank, Nx_box, Ny_box, Ny, nbx_box, nby_box, ncpu, m_tot, rho_tot, rho_loc, m_loc, rho_rec, m_rec);

  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    printf("Scattering done\n");

  
  /* printMatBorders(ncpu, rank, rho_loc, m_loc, m_rec, Nx_box, Ny_box, nbx_box, nby_box, Nx, Ny); */

  /* Send the borders for the cpus to communicate */
  sendBorders(rank, m_loc, rec_hborder, frame, Nx_box, Ny_box, sent_hborder);
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    printf("OK\n");
  sendBorders(rank, rho_loc, rec_hborder, frame, Nx_box, Ny_box, sent_hborder);

  
  if(rank==0)
    printf("Borders sent\n");

  /* printMatBorders(ncpu, rank, rho_loc, m_loc, m_rec, Nx_box, Ny_box, nbx_box, nby_box, Nx, Ny); */

  /* Print the initial field matrix */
  int t=0;
  gatherFields(ncpu, rank, rho_loc, m_loc, rho_tot, m_tot, m_rec0, rho_rec0, Nx_box, Ny_box, nbx_box, nby_box, Nx, Ny);

  if(rank==0)
    printf("Fields gathered\n");

  
  outputFields(rank, t, dt, Nx, Ny, f_rho, f_m, rho_tot, m_tot);

  
  if(rank==0)
    printf("Fields written\n");

  /* Initial position of the fluctuation's front */
  if(rank==0)
    {
      x0=0; // Absolute position of the fluctuation
      while(m_tot[x0][Ny2]*m_tot[x0+1][Ny2]>0)
	{
	  x0++;
	}
      box0=(Ny2/Ny_box)*nbx_box+floor(x0/Nx_box); // Box where to look for the front
      y0= Ny2%Ny_box;
    }

  MPI_Bcast(&x0, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&y0, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&box0, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  /* Temporal loop */
  for(t=0 ; t<Nt ; t++)
    {
      printf("%f\n", t*dt);

      /* updateFields(Nx_box, Ny_box, rho_loc, m_loc, rho_loc_old, m_loc_old, beta, gamma, c_adv, c_diff, dt); */

      updateFieldsNoTail(Nx_box, Ny_box, nbx_box, rho_loc, m_loc, rho_loc_old, m_loc_old, beta, gamma, c_adv, c_diff, dt, rank, box0, x0, ml, rhol, rank_x, shift, x0_shifted, box_x, x0_rel, x01, x01_rel, box_x1);

      /* Find position of the fluctuation */
      old_box0=box0;
      old_x0=x0;

      if(rank==box0)
	{
	  /* printf("x0_rel = %d\n", x0%Nx_box); */
	  if(m_loc[x0%Nx_box+1][y0]<0)
	    {
	      /* printf("OK 1\n"); */
	      if(x0%Nx_box==0)
		{
		  box0=frame[4*box0];
		}

	      x0--;

	      if(x0<0)
		x0=Nx-1;
	    }
	}

      MPI_Bcast(&x0, 1, MPI_INT, old_box0, MPI_COMM_WORLD);
      MPI_Bcast(&box0, 1, MPI_INT, old_box0, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);

      if(x0!=old_x0)
	{

	  rank_x = rank%nbx_box;
	    
	  x0_shifted = ((x0+Nx-shift)%(Nx));
	  box_x = (int) floor(x0_shifted/Nx_box);
	  x0_rel=x0_shifted%Nx_box;
  
	  x01 = ((x0+Nx-shift-1)%(Nx)); // Absolute position to the left of the fluctuation front
	  box_x1 = (int) floor(x01/Nx_box);
	  x01_rel = x01%Nx_box; // Relative position on the left site of the fluctuation front

	  if(rank_x==box_x1)
	    {

	      for(y=1 ; y<Ny_box+1 ; y++)
		{
		  rho_loc[x01_rel+1][y]=rhol;
		  m_loc[x01_rel+1][y]=ml;
		}
	    }

	}
      
      sendBorders(rank, m_loc, rec_hborder, frame, Nx_box, Ny_box, sent_hborder);
      sendBorders(rank, rho_loc, rec_hborder, frame, Nx_box, Ny_box, sent_hborder);

      /* Output snapshots */
      if(t>=clap)
	{
	  clap+=Nout;
	  gatherFields(ncpu, rank, rho_loc, m_loc, rho_tot, m_tot, m_rec0, rho_rec0, Nx_box, Ny_box, nbx_box, nby_box, Nx, Ny);
  
	  outputFields(rank, t, dt, Nx, Ny, f_rho, f_m, rho_tot, m_tot);

	  if(rank==0)
	    printf("%f\n", t*dt);
	  
	}
      
    }

  // CLOSE FILES
  fclose(f_input);
  
  // FREE ALLOC
  if(rank==0)
    {
      free_dmatrix(rho_tot);
      free_dmatrix(m_tot);
    }

  free_dmatrix(m_rec);
  free_dmatrix(rho_rec);

  free_dmatrix(m_rec0);
  free_dmatrix(rho_rec0);

  free_dmatrix(rho_loc);
  free_dmatrix(m_loc);
  free_dmatrix(rho_loc_old);
  free_dmatrix(m_loc_old);
  
  free(frame);
  free(rec_hborder);
  
  /* MPI_Type_free(&sent_hborder); */
  
  // CLOSE MPI
  MPI_Finalize();
  return(0);
}
