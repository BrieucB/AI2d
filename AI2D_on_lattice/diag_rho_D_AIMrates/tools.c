#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <float.h>
#include <omp.h>

#include "pcg_basic.h"
#include "tools.h"
#include "exit_if.h"
#include "erreurs.h"
#include "allocation.h"

double maxdouble32;
void constant(void)
{
  int i;
  
  maxdouble32=1.;
  
  for (i=0 ; i<32 ; i++)
    maxdouble32 *= 2.;
    
  return;
}


/* random double in [0;1] */
double genrand32_real1(pcg32_random_t *rng)
{
  return (((double) pcg32_random_r(rng))/(maxdouble32-1.)) ;
}

/* random double in [0;1[ */
double genrand32_real2(pcg32_random_t *rng)
{
  return (((double) pcg32_random_r(rng))/maxdouble32) ;
}

/* Compute fields */
double computeLocalQuantities(int lx, int ly, int N, int *x, int *y, int *s, int **local_m, int **local_rho)
{

  int tot_mag=0;
		
  /* RESET MATRICES */
  for(int xi=0 ; xi<lx ; xi++)
    {
      for(int yi=0 ; yi<ly ; yi++)
	{
	  local_m[yi][xi]=0;
	  local_rho[yi][xi]=0;
	}
    }

  
  /* FILL MATRICES */
  #pragma omp parallel for default(shared) reduction(+: tot_mag)
  for(int i=0 ; i<N ; i++)
    {
      int xi=x[i];
      int yi=y[i];
      int si=s[i];
      local_m[yi][xi]+=si;
      local_rho[yi][xi]+=1;
      tot_mag+=si;
    }

  return (((double) tot_mag)/((double) N));
}

/* Dynamics update */
void updatePositions(int lx, int ly, int *x, int *y, int *s, int N, int **local_m, int **local_rho, double w0, double beta, double eps, double D, double dt, pcg32_random_t *rng)
{

#pragma omp parallel for default(shared)
  for(int i=0 ; i<N ; i++)
    {
      int xi=x[i];
      int yi=y[i];
      int si=s[i];

      int mi=local_m[yi][xi];
      int rhoi=local_rho[yi][xi];
	  
      /* Flip the particle */
      double W=w0*exp(-si*beta*mi/rhoi);
      double p=genrand32_real2(rng);
      
      if(p<W*dt)
	{
	  s[i]=-si;
	}

      /* Move the particle to the right */
      else if(p<(W+D*(1.+s[i]*eps))*dt) 
	{
	  int x1=xi+1;
	  /* PBC in 1D */
	  if(x1>=lx)
	    x1-=lx;
	  x[i]=x1;
	}      

      /* Move the particle to the left */
      else if(p<(W+2.*D)*dt) 
	{
	  int x1=xi-1;
	  /* PBC in 1D */
	  if(x1<0)
	    x1+=lx;
	  x[i]=x1;
	}

      /* Diffuse up */
      else if(p<(W+3.*D)*dt)
	{
	  int y1 = yi+1;
	  /* PBC in 1D */
	  if(y1>=ly)
	    {
	      y1-=ly;
	    }
	  y[i]=y1;
	}

      /* Diffuse down */
      else if(p<(W+4.*D)*dt)
	{
	  int y1 = yi-1;
	  /* PBC in 1D */
	  if(y1<0)
	    {
	      y1+=ly;
	    }
	  y[i]=y1;
	}
    }
  

  return;
}
