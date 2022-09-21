#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <omp.h>
#include <time.h>

#include "allocation.h"
#include "exit_if.h"
#include "erreurs.h"

void solve_HD_out(int lx, int ly, int tmax, double dt, double ds, double beta, double D, double v, double rhol, double gamma, double tgap, double rhof)
{
  double Nt= floor(tmax/dt);

  int Nx= (int) floor(lx/ds) +2; // PADDING FOR PBC
  int Nx2= (int) floor(Nx/2);

  int Ny= (int) floor(ly/ds) +2;
  int Ny2= (int) floor(Ny/2);

  int Nout=(int) floor(tgap/dt);
  int clap = Nout;

  double c_diff=D*(dt/(ds*ds));
  double c_adv=(v/2.)*(dt/ds);

  //  clock_t start, end;
  double start, end;
  
  /* Determine ml = rhol*tanh(beta*phil) */
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


  int Nxfront=(int) floor(0.1*Nx);
  
  double r0=5;
  for(int x=Nxfront-((int) floor(r0/ds))+1 ; x<Nxfront+((int) floor(r0/ds)); x++)
    {
      double dsx=((double) x- (double) Nxfront);

      for(int y=Ny2-((int) floor(r0/ds))+1 ; y<Ny2+((int) floor(r0/ds)) ; y++)
	{
	  double dsy=((double) y - (double) Ny2);
	  if(ds*ds*(dsx*dsx+dsy*dsy)<=r0*r0)
	    {
	      rho[x][y] += (rhof/0.148)*exp(-1./(1.-ds*ds*(dsx*dsx+dsy*dsy)/(r0*r0)));
	      m[x][y] -= (rhof/0.148)*exp(-1./(1.-ds*ds*(dsx*dsx+dsy*dsy)/(r0*r0))); // ADD A BLOB OF density rhof : change normalization constant in 1d
	    }
	}
    }


  /* for(int x=Nx2-((int) floor(r0/ds))+1 ; x<Nx2+((int) floor(r0/ds)); x++) */
  /*   { */
  /*     double dsx=((double) x- (double) Nx2); */

  /*     for(int y=Ny2-((int) floor(r0/ds))+1 ; y<Ny2+((int) floor(r0/ds)) ; y++) */
  /* 	{ */
  /* 	  double dsy=((double) y - (double) Ny2); */
  /* 	  if(ds*ds*(dsx*dsx+dsy*dsy)<=r0*r0) */
  /* 	    { */
  /* 	      rho[x][y] += (rhof/0.148)*exp(-1./(1.-ds*ds*(dsx*dsx+dsy*dsy)/(r0*r0))); */
  /* 	      m[x][y] -= (rhof/0.148)*exp(-1./(1.-ds*ds*(dsx*dsx+dsy*dsy)/(r0*r0))); // ADD A BLOB OF density rhof : change normalization constant in 1d */
  /* 	    } */
  /* 	} */
  /*   } */
  
  FILE *f_rho, *f_m;

  /* f_rho=fopen("f_rho_t0.dat", "w"); */
  /* f_m=fopen("f_m_t0.dat", "w"); */


  /* for(int x=1 ; x<Nx2 ; x++) */
  /*   { */

  /*     for(int y=1 ; y<Ny-1 ; y++) */
  /* 	{ */

  /* 	  fprintf(f_rho, "%f ", rho[x][y]); */
  /* 	  fprintf(f_m, "%f ", m[x][y]); */

  /* 	} */
      
  /*     fprintf(f_rho, "\n"); */
  /*     fprintf(f_m, "\n"); */

  /*   } */
  
  /* fprintf(f_rho, "\n"); */
  /* fprintf(f_m, "\n"); */


  /* int xmin = Nx2-100; */
  /* int xmax = Nx2+100; */
  /* int ymin = Ny2-100; */
  /* int ymax = Ny2+100; */

  double t=0;
  int x0=0;

  
  int t1=0;

  start=omp_get_wtime();//clock();
  while(t<Nt)
    {
      //      printf("%.2f\n", t*dt);

      if((int) floor(t*dt)>=t1)
	{
	  end=omp_get_wtime();//clock();
	  double time_taken = ((double) end-start);
	  printf("%d %f SEC\n", (int) floor(t*dt), time_taken);
	  t1+=10;
	}
      
      /* for(int x=0; x<Nx ; x++) */
      /* 	{ */
  
      /* 	  memcpy(old_rho[x], rho[x], Ny*sizeof(double)); */
      /* 	  memcpy(old_m[x], m[x], Ny*sizeof(double)); */
      /* 	} */

      int x=0;
      while(m[x][Ny2]*m[x+1][Ny2]>0)
	{
	  x++;
	}
      int xmin = x;

      //      printf("xmin  = %d\n", xmin);

      /* Roll the lattice along the x-axis to set the droplet front on the left boundary */


      //printf("xmax  = %d\n", Nx+Nxfront-xmin-1);
      
#pragma omp parallel for default(shared)      
      for(int x=0 ; x<Nx ; x++)
	{
	  /* printf("u = %d ", x); */
	  /* printf("x = %d \n", x-Nxfront+xmin); */

	  for(int y=0 ; y<Ny ; y++)
	    {
	      old_m[x][y] = m[(x-Nxfront+xmin+Nx)%Nx][y];
	      old_rho[x][y] = rho[(x-Nxfront+xmin+Nx)%Nx][y];
	    }

	}

      x0+=-Nxfront+xmin;
      
      /* Fill the blank with liquid phase */
      for(int x=Nx-10 ; x<Nx ; x++)
	{
	  for(int y=0 ; y<Ny ; y++)
	    {
	      old_m[x][y]=ml;
	      old_rho[x][y]=rhol;
	    }
	}

      for(int x=0 ; x<10 ; x++)
	{
	  for(int y=0 ; y<Ny ; y++)
	    {
	      old_m[x][y]=ml;
	      old_rho[x][y]=rhol;
	    }
	} 


      /* Roll the lattice along the x-axis to center the droplet
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

#pragma omp parallel for default(shared)
      /* Update density & magnetization */
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
			  - 4*old_rho[x][y]);

	      m[x][y] = old_m[x][y]
		- c_adv*(old_rho[x+1][y]-old_rho[x-1][y])
		+ c_diff*(old_m[x-1][y]
			  + old_m[x+1][y]
			  + old_m[x][y-1]
			  + old_m[x][y+1]
			  - 4*old_m[x][y])
		+ 2.*dt*gamma*(old_rho[x][y]*sinh(betaphi)
			       - old_m[x][y]*cosh(betaphi));
	    }
	}

      /* PADDING */
      // Top & bottom rows
      for(int x=0 ; x<Nx ; x++)
	{

	  //	  printf("%d \n", x);
	  rho[x][0] = rho[x][Ny-2];
	  m[x][0] = m[x][Ny-2];

	  //	  printf("rho[%d][1] = %lg\n", x, rho[x][1]);
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
	  /* printf("%f\n", t*dt); */

	  /* /\* Find horizontal range spanned by the droplet *\/ */
	  /* int x=0; */
	  /* while(m[x][Ny2]*m[x+1][Ny2]>0) */
	  /*   { */
	  /*     x++; */
	  /*   } */
	  /* int xmin = x; */
	  /* //	  printf("xmin = %d\n", xmin); */

	  /* x+=10; */
	  /* while(m[x][Ny2]<-1e-1) */
	  /*   { */
	  /*     x++; */
	  /*   } */
	  /* int xmax = x; */
	  /* //	  printf("xmax = %d\n", xmax);	   */

	  /* /\* Find vertical range spanned by the droplet *\/	   */
	  /* int Lmax=0; */

	  /* for(int x=xmin ; x<xmax ; x++) */
	  /*   { */
	  /*     int sum=0; */
	  /*     for(int y=1 ; y<Ny ; y++) */
	  /* 	{ */
	  /* 	  if(m[x][y]<0) */
	  /* 	    sum+=1; */
	  /* 	} */

	  /*     if(sum>Lmax) */
	  /* 	Lmax=sum; */
	  /*   } */

	  /* int ymin = Ny2 - ((int) floor(Lmax/2)); */
	  /* int ymax = Ny2 + ((int) floor(Lmax/2)); */
	  /* //	  printf("Lmax = %d\nymin = %d\nymax = %d\n", Lmax, ymin, ymax);	  	   */

	  /* Create output files */
	  char f_name_m[200];
	  int pos_m = 0;
	  pos_m += sprintf(&f_name_m[pos_m], "f_m_t%d.dat", (int) floor(t*dt));
	  f_m=fopen(f_name_m, "w");

	  char f_name_rho[200];
	  int pos_rho = 0;
	  pos_rho += sprintf(&f_name_rho[pos_rho], "f_rho_t%d.dat", (int) floor(t*dt));
	  f_rho=fopen(f_name_rho, "w");

	  
	  for(int x=1 ; x<Nx-1 ; x++) //for(int x=(xmin-100)%Nx ; x<(xmax+100)%Nx ; x++)
	    {
	      fprintf(f_rho, "%d ", x+x0);
	      fprintf(f_m, "%d ", x+x0);

	      for(int y=1 ; y<Ny-1 ; y++) //for(int y=(ymin-100)%Ny ; y<(ymax+100)%Ny ; y++)
		{
		  /* char *rho_xy=(char*)malloc(50 * sizeof(char)); */
		  /* sprintf(rho_xy, "%f", rho[x][y]); */
		  /* printf("%s ", rho_xy); */

		  //printf("%f ", rho_xy);
		  fprintf(f_rho, "%f ", rho[x][y]);
		  fprintf(f_m, "%f ", m[x][y]);

		}
      
	      fprintf(f_rho, "\n");
	      fprintf(f_m, "\n");

	    }
  
	  fprintf(f_rho, "\n");
	  fprintf(f_m, "\n");

	  fclose(f_rho);
	  fclose(f_m);

	}
	
      t++;
    }
    
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

  solve_HD_out(lx, ly, tmax, dt, ds, beta, D, v, rhol, gamma, tgap, rhof);

  fclose(f_input);

  return(0);
}
