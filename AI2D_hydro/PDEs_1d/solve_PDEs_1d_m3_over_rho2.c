#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <omp.h>

void solve_HD_out(int L, int tmax, double dt, double dx, double beta, double r, double D, double v, double rhol, double tgap, double rhoc)
{
  int Nt= (int) floor(tmax/dt);
  int Nx= (int) floor(L/dx);
  int Nx2= (int) floor(Nx/2);

  int Nout=(int) floor(tgap/dt);
  int clap = 0;
  
  double c_diff=D*(dt/(dx*dx));
  double c_adv=0.5*v*(dt/dx);
  double c_m=2*(beta-1);
  double c_m3=dt*beta*beta*(1-beta/3);
  
  double ml=rhol*sqrt(2*(beta-1)/(beta*beta*(1-beta/3)));
  double mc=rhoc;

  double *rho = (double *)malloc(Nx*sizeof(double));
  double *m = (double *)malloc(Nx*sizeof(double));
  double *old_rho = (double *)malloc(Nx*sizeof(double));
  double *old_m = (double *)malloc(Nx*sizeof(double));

  for(int x=0 ; x<Nx ; x++)
    {
      rho[x]=rhol;
      m[x]=ml;
    }

  double r0=5;
  for(int x=Nx2-((int) floor(r0/dx))+1 ; x<Nx2+((int) floor(r0/dx)) ; x++)
    {
      double dx1=((double) x- (double) Nx2);
      rho[x] += rhoc*exp(-1./(1-dx*dx*dx1*dx1/(r0*r0)));
      m[x] -= 3.*mc*exp(-1./(1-dx*dx*dx1*dx1/(r0*r0))); // ADD A BLOB OF density rho
    }
  
  FILE *f_rho, *f_m;

  f_rho=fopen("f_rho.dat", "w");
  f_m=fopen("f_m.dat", "w");

  for(int x=0 ; x<Nx-1 ; x++)
    {

      fprintf(f_rho, "%d ", x);
      fprintf(f_m, "%d ", x);
	  
    }
      
  fprintf(f_rho, "%d\n", Nx-1);
  fprintf(f_m, "%d\n", Nx-1);


  for(int x=0 ; x<Nx-1 ; x++)
    {

      fprintf(f_rho, "%f ", rho[x]);
      fprintf(f_m, "%f ", m[x]);
	  
    }
  fprintf(f_rho, "%f\n", rho[Nx-1]);
  fprintf(f_m, "%f\n", m[Nx-1]);

  int t=0;
  while(t<Nt)
    {
      
      memcpy(old_rho, rho, Nx*sizeof(double));
      memcpy(old_m, m, Nx*sizeof(double));

      rho[0]=old_rho[0] - c_adv*(old_m[1]-old_m[Nx-1]) + c_diff*(old_rho[Nx-1]+old_rho[1]-2*old_rho[0]);
      m[0]=old_m[0] - c_adv*(old_rho[1]-old_rho[Nx-1]) + c_diff*(old_m[Nx-1]+old_m[1]-2*old_m[0]) + (c_m - 2.*r/old_rho[0])*dt*old_m[0] - c_m3*old_m[0]*old_m[0]*old_m[0]/(old_rho[0]*old_rho[0]);
      
      rho[Nx-1]=old_rho[Nx-1] - c_adv*(old_m[0]-old_m[Nx-2]) + c_diff*(old_rho[Nx-2]+old_rho[0]-2*old_rho[Nx-1]);
      m[Nx-1]=old_m[Nx-1] - c_adv*(old_rho[0]-old_rho[Nx-2]) + c_diff*(old_m[Nx-2]+old_m[0]-2*old_m[Nx-1]) + (c_m - 2.*r/old_rho[Nx-1])*dt*old_m[Nx-1] - c_m3*old_m[Nx-1]*old_m[Nx-1]*old_m[Nx-1]/(old_rho[Nx-1]*old_rho[Nx-1]);
      
      for(int x=1 ; x<Nx-1 ; x++)
	{
	  rho[x]=old_rho[x] - c_adv*(old_m[x+1]-old_m[x-1]) + c_diff*(old_rho[x-1]+old_rho[x+1]-2*old_rho[x]);
	  m[x]=old_m[x] - c_adv*(old_rho[x+1]-old_rho[x-1]) + c_diff*(old_m[x-1]+old_m[x+1]-2*old_m[x]) + (c_m - 2.*r/old_rho[x])*dt*old_m[x] - c_m3*old_m[x]*old_m[x]*old_m[x]/(old_rho[x]*old_rho[x]);

	  /* if(rho[x]<DBL_EPSILON) */
	  /*   { */
	  /*     rho[x]=0.; */
	  /*     m[x]=0.; */
	  /*   } */
						 
	}

      if(t>=clap)
	{
	  clap+=Nout;
	  
	  for(int x=0 ; x<Nx-1 ; x++)
	    {

	      fprintf(f_rho, "%f ", rho[x]);
	      fprintf(f_m, "%f ", m[x]);
      
	    }
	  fprintf(f_rho, "%f\n", rho[Nx-1]);
	  fprintf(f_m, "%f\n", m[Nx-1]);
	}
	
      t++;
    }
  
  fclose(f_rho);
  fclose(f_m);
    
  free(rho);
  free(m);
  free(old_rho);
  free(old_m);
}

int main(void)
{

  FILE *f_input;
  f_input=fopen("f_input.dat", "r");

  int L;
  int tmax;
  double tgap;

  
  double dt;
  double dx;
  double beta;
  double D;
  double v;
  double rhol;
  double r;

  double rhoc;

  if(fscanf(f_input, "tgap = %lg tmax = %d dt = %lg lx = %d dx = %lg rhol = %lg beta = %lg v = %lg D = %lg r = %lg rhoc = %lg", &tgap, &tmax, &dt, &L, &dx, &rhol, &beta, &v, &D, &r, &rhoc)!=11){exit(1);};

  solve_HD_out(L, tmax, dt, dx, beta, r, D, v, rhol, tgap, rhoc);

  return(0);
}
