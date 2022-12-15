#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <omp.h>

#include "allocation.h"
#include "exit_if.h"
#include "erreurs.h"

void solve_HD_out(int lx, int ly, int tmax, double dt, double ds, double beta, double D, double v, double gamma, double rhol, double tgap, double rhof)
{
  int Nt= (int) floor(tmax/dt);

  int Nx= (int) floor(lx/ds) +2; // PADDING FOR PBC
  int Nx2= (int) floor(Nx/2);

  int Ny= (int) floor(ly/ds) +2;
  int Ny2= (int) floor(Ny/2);

  int Nout=(int) floor(tgap/dt);
  int clap = 0;

  double c_diff=D*(dt/(ds*ds));
  double c_adv=(v/2.)*(dt/ds);
  
  double phil = 0.5;
  for(int p=0 ; p<100 ; p++)
    {
      phil=tanh(beta*phil);
    }
  double ml=rhol*phil;
    
  double **rho = dmatrix((long) Nx, (long) Ny);
  double **m = dmatrix((long) Nx, (long) Ny);

  double **old_rho = dmatrix((long) Nx, (long) Ny);
  double **old_m = dmatrix((long) Nx, (long) Ny);

  for(int x=0 ; x<Nx ; x++)
    {
      for(int y=0 ; y<Ny ; y++)
	{

	  rho[x][y]=rhol;
	  m[x][y]=ml;


	  old_rho[x][y]=rhol;
	  old_m[x][y]=ml;

	}
    }

  printf("Nx2 = %d, rhol = %f, ml = %f\n", Nx2, rhol, ml);
  
  double r0=2;
  for(int x=Nx2-((int) floor(r0/ds))+1 ; x<Nx2+((int) floor(r0/ds)); x++)
    {
      double dsx=((double) x- (double) Nx2);

      for(int y=Ny2-((int) floor(r0/ds))+1 ; y<Ny2+((int) floor(r0/ds)) ; y++)
	{
	  double dsy=((double) y - (double) Ny2);
	  /* printf("x = %d, y = %d\n", x, y); */
	  if(ds*ds*(dsx*dsx+dsy*dsy)<=r0*r0)
	    {
	      rho[x][y] += (rhof/0.148)*exp(-1./(1.-ds*ds*(dsx*dsx+dsy*dsy)/(r0*r0)));
	      m[x][y] -= (rhof/0.148)*exp(-1./(1.-ds*ds*(dsx*dsx+dsy*dsy)/(r0*r0))); // ADD A BLOB OF density rhof : change normalization constant in 1d
/* 	      printf("val = %f\n", (rhof/0.148)*exp(-1./(1.-ds*ds*(dsx*dsx+dsy*dsy)/(r0*r0)))); */
/* 	      printf("val_m = %f\n\n",ml- (rhof/0.148)*exp(-1./(1.-ds*ds*(dsx*dsx+dsy*dsy\ */
/* )/(r0*r0)))); */
	    }
	  
	}
      /* printf("\n\n"); */
    }

  
  FILE *f_rho, *f_m;

  f_rho=fopen("f_rho.dat", "w");
  f_m=fopen("f_m.dat", "w");


  for(int x=1 ; x<Nx-1 ; x++)
    {

      for(int y=1 ; y<Ny-1 ; y++)
	{

	  fprintf(f_rho, "%f ", rho[x][y]);
	  fprintf(f_m, "%f ", m[x][y]);

	}
      
      fprintf(f_rho, "\n");
      fprintf(f_m, "\n");

    }
  
  fprintf(f_rho, "\n");
  fprintf(f_m, "\n");

      
  int t=0;
  clap+=Nout;
  while(t<Nt)
    {
      //#pragma omp parallel for default(shared)      
      /* for(int x=0; x<Nx ; x++) */
      /* 	{ */
  
      /* 	  memcpy(old_rho[x], rho[x], Ny*sizeof(double)); */
      /* 	  memcpy(old_m[x], m[x], Ny*sizeof(double)); */
      /* 	} */

      /* memcpy(old_rho[0], rho[0], (Nx*Ny)*sizeof(double)); */
      /* memcpy(old_m[0], m[0], (Nx*Ny)*sizeof(double)); */
      /* printf("dt = %lg beta = %lg gamma = %lg c_adv = %lg c_diff = %lg\n", dt, beta, gamma, c_adv, c_diff); */

      for(int x=0 ; x<Nx ; x++)
	{
	  for(int y=0 ; y<Ny ; y++)
	    {
	      old_m[x][y]=m[x][y];
	      old_rho[x][y]=rho[x][y];
	    }
	}
      
#pragma omp parallel for default(shared)
      for(int x=1 ; x<Nx-1 ; x++)
	{
	  for(int y=1 ; y<Ny-1 ; y++)
	    {
	      double betaphi = beta*old_m[x][y]/old_rho[x][y];
	      
	      rho[x][y] = old_rho[x][y]
		- c_adv*(old_m[x+1][y]-old_m[x-1][y]) // Advection                        
		+ c_diff*(old_rho[x-1][y] // Diffusion                                        
			  + old_rho[x+1][y]
			  + old_rho[x][y-1]
			  + old_rho[x][y+1]
			  - 4.*old_rho[x][y]);

	      m[x][y] = old_m[x][y]
		- c_adv*(old_rho[x+1][y]-old_rho[x-1][y])
		+ c_diff*(old_m[x-1][y]
			  + old_m[x+1][y]
			  + old_m[x][y-1]
			  + old_m[x][y+1]
			  - 4.*old_m[x][y])
		+ 2.*dt*gamma*(old_rho[x][y]*sinh(betaphi)
			       - old_m[x][y]*cosh(betaphi));
	      /* printf("%d %d %f\n", x, y, m[x][y]); */
	    }
	}
      /* for(int x=1 ; x<Nx-1 ; x++) */
      /* 	{ */
      /* 	  for(int y=1 ; y<Ny-1 ; y++) */
      /* 	    { */

      /* 	      rho[x][y] = old_rho[x][y] - c_adv*(old_m[x+1][y]-old_m[x-1][y]) + c_diff*(old_rho[x-1][y] + old_rho[x+1][y]  + old_rho[x][y-1]  + old_rho[x][y+1]  - 4*old_rho[x][y]); */

      /* 	      m[x][y] = old_m[x][y] */
      /* 		- c_adv*(old_rho[x+1][y]-old_rho[x-1][y]) */
      /* 		+ c_diff*(old_m[x-1][y] */
      /* 			  + old_m[x+1][y] */
      /* 			  + old_m[x][y-1] */
      /* 			  + old_m[x][y+1] */
      /* 			  - 4*old_m[x][y]) */
      /* 		+ (c_m - 2.*r/old_rho[x][y])*dt*old_m[x][y] */
      /* 		- c_m3*old_m[x][y]*old_m[x][y]*old_m[x][y]/(old_rho[x][y]*old_rho[x][y]); */
	  
	  /*     /\* if(rho[x]<DBL_EPSILON) *\/ */
	  /*     /\*   { *\/ */
	  /*     /\*     rho[x]=0.; *\/ */
	  /*     /\*     m[x]=0.; *\/ */
	  /*     /\*   } *\/ */
						 
	/*     } */
	/* } */

      /* PADDING */
      // Top & bottom rows
      for(int x=0 ; x<Nx ; x++)
	{
	  rho[x][0] = rho[x][Ny-2];
	  m[x][0] = m[x][Ny-2];

	  rho[x][Ny-1] = rho[x][1];
	  m[x][Ny-1] = m[x][1];
	}
      
      // Left & right cols
      for(int y=1 ; y<Ny-1 ; y++)
	{
	  rho[0][y] = rho[Nx-2][y];
	  m[0][y] = m[Nx-2][y];

	  rho[Nx-1][y] = rho[1][y];
	  m[Nx-1][y] = m[1][y];
	}

      if(t>=clap)
	{
	  clap+=Nout;
	  printf("%d\n", t);
	  for(int x=1 ; x<Nx-1 ; x++)
	    {

	      for(int y=1 ; y<Ny-1 ; y++)
		{
		  /* char *rho_xy=(char*)malloc(50 * sizeof(char)); */
		  /* sprintf(rho_xy, "%f", rho[x][y]); */
		  /* printf("%s ", rho_xy); */

		  double rho_xy=rho[x][y];
		  if(isinf(rho_xy))
		    {
		      exit(0);
		    }
			  
		  else
		    {
		      //printf("%f ", rho_xy);
		      fprintf(f_rho, "%f ", rho[x][y]);
		      fprintf(f_m, "%f ", m[x][y]);
		    }

		}
      
	      fprintf(f_rho, "\n");
	      fprintf(f_m, "\n");

	    }
  
	  fprintf(f_rho, "\n");
	  fprintf(f_m, "\n");
	}
	
      t++;
    }
  
  fclose(f_rho);
  fclose(f_m);
    
  free_dmatrix(rho);
  free_dmatrix(m);

  free_dmatrix(old_rho);
  free_dmatrix(old_m);
}

int main(void)
{

  FILE *f_input;
  f_input=fopen("f_input.dat", "r");

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

  fscanf(f_input, "tgap = %lg tmax = %d dt = %lg lx = %d ly = %d ds = %lg rhol = %lg beta = %lg v = %lg D = %lg gamma = %lg rhof = %lg", &tgap, &tmax, &dt, &lx, &ly, &ds, &rhol, &beta, &v, &D, &gamma, &rhof);

  printf("tgap = %lg tmax = %d dt = %lg lx = %d ly = %d ds = %lg rhol = %lg beta = %lg v = %lg D = %lg gamma = %lg rhof = %lg", tgap, tmax, dt, lx, ly, ds, rhol, beta, v, D, gamma, rhof);

  
  solve_HD_out(lx, ly, tmax, dt, ds, beta, D, v, gamma, rhol, tgap, rhof);

  fclose(f_input);

  return(0);
}
