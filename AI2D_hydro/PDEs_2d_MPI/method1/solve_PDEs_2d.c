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
  
  // SPACE-TIME DISCRETIZATION
  double Nt= floor(tmax/dt);
  int Nx= (int) floor(lx/ds);
  //  int Nx2= (int) floor(Nx/2);
  int Ny= (int) floor(ly/ds);
  int Ny2= (int) floor(Ny/2);
  int Nout=(int) floor(tgap/dt);
  int clap = Nout;
  double c_diff=D*(dt/(ds*ds));
  double c_adv=(v/2.)*(dt/ds);

  // MPI BOXES
  int nbx_box = (int) ncpu/nby_box; // ncpu is a multiple of nby_box
  // => nby_box*nbx_box = ncpu
  
  int Nx_box = (int) Nx/nbx_box;
  int Ny_box = (int) Ny/nby_box;

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
  double **rho_rec0;

  /* Vector to receive the horizontal borders */
  double *rec_hborder =  (double *)malloc((Nx_box+2)*sizeof(double));

  int *frame = (int *)malloc(4*nbx_box*nby_box*sizeof(int));
  init_frame(nbx_box, nby_box, frame);
  
  /* if(rank==0) */
  /*   { */

  /*     for(int box=0 ; box<nbx_box*nby_box ; box++) */
  /* 	{ */
  /* 	  printf("Nei of box %d :\n", box); */
  /* 	  for(int nei=0 ; nei<4 ; nei++) */
  /* 	    { */
  /* 	      printf("%d ", frame[4*box+nei]); */
  /* 	    } */
  /* 	  printf("\n\n"); */
  /* 	} */
	
  /*     printf("\n"); */
  /*   } */
  
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
      m_rec = dmatrix((long) Nx_box+2, (long) Ny_box+2);
    }

  /* Initialize fields */
  initFields(rank, rho_tot, m_tot, m_rec, Nx, Ny, Ny2, Nx_box, Ny_box, ds, rhof, rhol,ml);

  /* Scatter fields among the cpus */
  scatterTotalMatrices(rank, Nx_box, Ny_box, Ny, nbx_box, nby_box, ncpu, m_tot, rho_tot, rho_loc, m_loc, rho_rec, m_rec);
  
  /* printMatBorders(ncpu, rank, rho_loc, m_loc, m_rec, Nx_box, Ny_box, nbx_box, nby_box, Nx, Ny); */

  /* Send the borders for the cpus to communicate */
  sendBorders(rank, m_loc, rec_hborder, frame, Nx_box, Ny_box);
  sendBorders(rank, rho_loc, rec_hborder, frame, Nx_box, Ny_box);

  /* printMatBorders(ncpu, rank, rho_loc, m_loc, m_rec, Nx_box, Ny_box, nbx_box, nby_box, Nx, Ny); */

  /* Print the initial field matrix */
  int t=0;
  gatherFields(ncpu, rank, rho_loc, m_loc, rho_tot, m_tot, m_rec0, Nx_box, Ny_box, nbx_box, nby_box, Nx, Ny);
  
  outputFields(rank, t, dt, Nx, Ny, f_rho, f_m, rho_tot, m_tot);
  
  /* Initial position of the fluctuation's front */
  int x0, y0, box0;
  if(rank==0)
    {
      x0=0; // Absolute position of the fluctuation
      while(m_tot[x0][Ny2]*m_tot[x0+1][Ny2]>0)
	{
	  x0++;
	}
      box0=(Ny2/Ny_box)*nbx_box+floor(x0/Nx_box); // Box where to look for the front
      y0= Ny2%Ny_box;
      printf("Box0 = %d, y0 = %d\n", box0, y0);
    }

  MPI_Bcast(&x0, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&y0, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&box0, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  int old_box0;
  /* int old_x0; */
  
  for(t=0 ; t<Nt ; t++)
    {
      /* if(rank==0) */
      /* 	printf("%d\n", t); */
      
      /* Find position of the fluctuation */
      old_box0=box0;
      /* old_x0=x0; */

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

	      if(x0==0)
		x0=Nx;
	    }
	}
      /* if(rank==0) */
      /* 	printf("OK\n"); */
      MPI_Bcast(&x0, 1, MPI_INT, old_box0, MPI_COMM_WORLD);
      MPI_Bcast(&box0, 1, MPI_INT, old_box0, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);

      /* if(x0!=old_x0) */
      /* 	{ */

      /* 	  int rank_x = rank%nbx_box; */
      /* 	  int shift=10; */
  
      /* 	  int x0_shifted = ((x0+Nx-shift)%(Nx)); */
      /* 	  int box_x = (int) floor(x0_shifted/Nx_box); */
      /* 	  int x0_rel=x0_shifted%Nx_box; */
  
      /* 	  int x01 = ((x0+Nx-shift-1)%(Nx)); // Absolute position to the left of the fluctuation front */
      /* 	  int box_x1 = (int) floor(x01/Nx_box); */
      /* 	  int x01_rel = x01%Nx_box; // Relative position on the left site of the fluctuation front */

      /* 	  if(rank_x==box_x1) */
      /* 	    { */

      /* 	      for(int y=1 ; y<Ny_box+1 ; y++) */
      /* 		{ */
      /* 		  rho_loc[x01_rel+1][y]=rhol; */
      /* 		  m_loc[x01_rel+1][y]=ml; */
      /* 		} */
      /* 	    } */

	        
      /* 	  sendBorders(rank, m_loc, rec_hborder, frame, Nx_box, Ny_box); */
      /* 	  sendBorders(rank, rho_loc, rec_hborder, frame, Nx_box, Ny_box); */


      /* 	} */
      
      /* if(rank==0) */
      /* 	printf("box0 = %d, x0 = %d\n", box0, x0); */
      
      /* updateFields(Nx_box, Ny_box, rho_loc, m_loc, rho_loc_old, m_loc_old, beta, gamma, c_adv, c_diff, dt); */

      updateFieldsNoTail2(Nx_box, Ny_box, nbx_box, rho_loc, m_loc, rho_loc_old, m_loc_old, beta, gamma, c_adv, c_diff, dt, rank, box0, x0, ml, rhol);
      
      sendBorders(rank, m_loc, rec_hborder, frame, Nx_box, Ny_box);
      sendBorders(rank, rho_loc, rec_hborder, frame, Nx_box, Ny_box);

      if(t>=clap)
	{
	  clap+=Nout;
	  gatherFields(ncpu, rank, rho_loc, m_loc, rho_tot, m_tot, m_rec0, Nx_box, Ny_box, nbx_box, nby_box, Nx, Ny);
  
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
      free_dmatrix(m_rec);
    }

  free_dmatrix(rho_loc);
  free_dmatrix(m_loc);
  free_dmatrix(rho_loc_old);
  free_dmatrix(m_loc_old);
  
  free(frame);

  // CLOSE MPI
  MPI_Finalize();
  return(0);
}




/* int x=0; */
/* while(m[x][Ny2]*m[x+1][Ny2]>0) */
/* 	{ */
/* 	  x++; */
/* 	} */
/* int xmin = x; */

//      printf("xmin  = %d\n", xmin);

/* Roll the lattice along the x-axis to set the droplet front on the left boundary */


//printf("xmax  = %d\n", Nx+Nxfront-xmin-1);
      
/* for(int x=0 ; x<Nx ; x++) */
/* 	{ */
/* 	  /\* printf("u = %d ", x); *\/ */
/* 	  /\* printf("x = %d \n", x-Nxfront+xmin); *\/ */

/* 	  for(int y=0 ; y<Ny ; y++) */
/* 	    { */
/* 	      old_m[x][y] = m[(x-Nxfront+xmin+Nx)%Nx][y]; */
/* 	      old_rho[x][y] = rho[(x-Nxfront+xmin+Nx)%Nx][y]; */
/* 	    } */

/* 	} */

/* x0+=-Nxfront+xmin; */
      
/* /\* Fill the blank with liquid phase *\/ */
/* for(int x=Nx-10 ; x<Nx ; x++) */
/* 	{ */
/* 	  for(int y=0 ; y<Ny ; y++) */
/* 	    { */
/* 	      old_m[x][y]=ml; */
/* 	      old_rho[x][y]=rhol; */
/* 	    } */
/* 	} */

/* for(int x=0 ; x<10 ; x++) */
/* 	{ */
/* 	  for(int y=0 ; y<Ny ; y++) */
/* 	    { */
/* 	      old_m[x][y]=ml; */
/* 	      old_rho[x][y]=rhol; */
/* 	    } */
/* 	}  */


/* Roll the lattice along the x-axis to center the droplet */
/* for(int x=0 ; x<Nx2+xmin ; x++) */
/* 	{ */
/* 	  for(int y=1 ; y<Ny-1 ; y++) */
/* 	    { */
/* 	      old_m[x+Nx2-xmin][y]=m[x][y]; */
/* 	      old_rho[x+Nx2-xmin][y]=rho[x][y]; */
/* 	    } */
/* 	} */

/* /\* Fill the blank with liquid phase *\/ */
/* for(int x=0 ; x<3 ; x++) */
/* 	{ */
/* 	  for(int y=0 ; y<Ny ; y++) */
/* 	    { */
/* 	      old_m[x][y]=ml; */
/* 	      old_rho[x][y]=rhol; */
/* 	    } */
/* 	} */
