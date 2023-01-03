#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <omp.h>

void solve_HD_out(int lx, int tmax, double dt, double ds, double beta, double D, double v, double rhol, double gamma, double tgap, double rhof)
{
  int Nt= (int) floor(tmax/dt);

  int Nx= (int) floor(lx/ds) +2; // PADDING FOR PBC
  int Nx2= (int) floor(Nx/2);

  int Nout=(int) floor(tgap/dt);
  int clap = 0;//Nout;

  double c_diff=D*(dt/(ds*ds));
  double c_adv=(v/2.)*(dt/ds);

  /* Determine ml = rhol*tanh(beta*phil) */
  double phil = 0.5;

  for(int p=0 ; p<100 ; p++)
    {
      phil=tanh(beta*phil);
    }

  double ml=rhol*phil;
  /* double ml=rhol*sqrt(2*(beta-1)/(beta*beta*(1-beta/3))); */
    
  double *rho = (double *)malloc(Nx*sizeof(double));
  double *m = (double *)malloc(Nx*sizeof(double));

  double *old_rho = (double *)malloc(Nx*sizeof(double));
  double *old_m = (double *)malloc(Nx*sizeof(double));
  
  for(int x=0 ; x<Nx ; x++)
    {
      rho[x]=rhol;
      m[x]=ml;

      old_rho[x]=rhol;
      old_m[x]=ml;
    }


  /* int Nxfront=(int) floor(0.1*Nx); */
  int Nxfront=(int) floor(0.5*Nx);
  
  double r0=5;
  for(int x=Nxfront-((int) floor(r0/ds))+1 ; x<Nxfront+((int) floor(r0/ds)); x++)
    {
      double dsx=((double) x- (double) Nxfront);
      if(ds*ds*(dsx*dsx)<=r0*r0)
	{
	  rho[x] += (rhof/0.148)*exp(-1./(1.-ds*ds*(dsx*dsx)/(r0*r0)));
	  m[x] -= (rhof/0.148)*exp(-1./(1.-ds*ds*(dsx*dsx)/(r0*r0))); // ADD A BLOB OF density rhof : change normalization constant in 1d
	}
    }

  FILE *f_rho, *f_m;
  FILE *f_out;

  f_out=fopen("f_out.dat", "w");
  
  int t=0;
  int x0=0;

  int t1=0;

  /* int cpt_t=0; */


  while(t<Nt)
    {
      /* printf("%.2f\n", t*dt); */

      /* int x=0; */
      /* while(m[x]*m[x+1]>0) */
      /* 	{ */
      /* 	  x++; */
      /* 	} */
      /* int xmin = x; */

      /* Roll the lattice along the x-axis to set the droplet front on the left boundary */
      
/* #pragma omp parallel for default(shared)       */
/*       for(int x=0 ; x<Nx ; x++) */
/* 	{ */
/* 	  old_m[x] = m[(x-Nxfront+xmin+Nx)%Nx]; */
/* 	  old_rho[x] = rho[(x-Nxfront+xmin+Nx)%Nx]; */
/* 	} */

/*       x0+=-Nxfront+xmin; */
      
/*       /\* Fill the blank with liquid phase *\/ */
/*       for(int x=Nx-10 ; x<Nx ; x++) */
/* 	{ */
/* 	  old_m[x]=ml; */
/* 	  old_rho[x]=rhol; */
/* 	} */

/*       for(int x=0 ; x<10 ; x++) */
/* 	{ */
/* 	  old_m[x]=ml; */
/* 	  old_rho[x]=rhol; */
/* 	}  */

            
#pragma omp parallel for default(shared)      
      for(int x=0 ; x<Nx ; x++)
	{
	  old_m[x] = m[x];
	  old_rho[x] = rho[x];
	}

      
#pragma omp parallel for default(shared)
      /* Update density & magnetization */
      for(int x=1 ; x<Nx-1 ; x++)
	{
	  double betaphi = beta*old_m[x]/old_rho[x];
	      
	  rho[x] = old_rho[x] 
	    - c_adv*(old_m[x+1]-old_m[x-1]) // Advection 
	    + c_diff*(old_rho[x-1] // Diffusion
		      + old_rho[x+1]
		      - 2*old_rho[x])
	    ;

	  m[x] = old_m[x]
	    - c_adv*(old_rho[x+1]-old_rho[x-1])
	    + c_diff*(old_m[x-1]
		      + old_m[x+1]
		      - 2*old_m[x])
	    /* - dt*old_m[x]*(old_m[x+1]-old_m[x-1])/2. */
	    + 2.*dt*gamma*(old_rho[x]*sinh(betaphi)
			   - old_m[x]*cosh(betaphi));
	  
	}

      /* PADDING */
      rho[0] = rho[Nx-2];
      m[0] = m[Nx-2];

      rho[Nx-1] = rho[1];
      m[Nx-1] = m[1];

      /* FRONT POSITION, BACK POSITION, FRONT DENSITY */

      if(t>=t1)
	{	

	  int frontpos=0;
	  int backpos=0;
	  int cometpos=0;
	  double frontdens=0;
	  double frontmag=0;
  
	  for(int x=1 ; x<Nx-2 ; x++)
	    {
	      if((m[x]<frontmag))
		frontmag=m[x];
	      
	    }
	  //	  printf("%f\n", frontmag);
	  
	  for(int x=1 ; x<Nx-2 ; x++)
	    {
	      if((m[x]>=0)&&(m[x+1]<0))
		frontpos=x;

	      if((m[x]<frontmag/2.)&&(m[x+1]>=frontmag/2.))
		backpos=x;

	      if((m[x]<0.5)&&(m[x+1]>=0.5))
		cometpos=x;
	  
	      if(rho[x]>frontdens)
		frontdens=rho[x];
	    }
	  t1+=(int) floor(10/dt);

	  fprintf(f_out, "%f %d %d %d %f\n", t*dt, frontpos, backpos, cometpos, frontdens);
	}

      
      
      if(t>=clap)
	{
	  clap+=Nout;
	  /* printf("%f\n", t*dt); */

	  /* Create output files */
	  char f_name_m[200];
	  int pos_m = 0;
	  pos_m += sprintf(&f_name_m[pos_m], "f_m_t%d.dat", (int) floor(t*dt));
	  /* pos_m += sprintf(&f_name_m[pos_m], "f_m_t%d.dat", cpt_t); */
	  f_m=fopen(f_name_m, "w");

	  char f_name_rho[200];
	  int pos_rho = 0;
	  pos_rho += sprintf(&f_name_rho[pos_rho], "f_rho_t%d.dat", (int) floor(t*dt));
	  /* pos_rho += sprintf(&f_name_rho[pos_rho], "f_rho_t%d.dat", cpt_t); */
	  f_rho=fopen(f_name_rho, "w");

	  /* cpt_t++; */
	  
	  for(int x=1 ; x<Nx-1 ; x++) //for(int x=(xmin-100)%Nx ; x<(xmax+100)%Nx ; x++)
	    {
	      fprintf(f_rho, "%d ", x+x0);
	      fprintf(f_m, "%d ", x+x0);

	      fprintf(f_rho, "%f\n", rho[x]);
	      fprintf(f_m, "%f\n", m[x]);

	    }
        
	  fclose(f_rho);
	  fclose(f_m);

	}
	
      t++;
    }
    
  free(rho);
  free(m);

  free(old_rho);
  free(old_m);
}

int main(void)
{

  FILE *f_input;
  f_input=fopen("f_input.dat", "r");

  int lx;
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

  fscanf(f_input, "tgap = %lg tmax = %d dt = %lg lx = %d ds = %lg rhol = %lg beta = %lg v = %lg D = %lg gamma = %lg rhof = %lg", &tgap, &tmax, &dt, &lx, &ds, &rhol, &beta, &v, &D, &gamma, &rhof);

  solve_HD_out(lx, tmax, dt, ds, beta, D, v, rhol, gamma, tgap, rhof);

  fclose(f_input);

  return(0);
}
