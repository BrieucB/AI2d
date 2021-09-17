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

int WhichBox(int ly, double x, double y) 
{
  int i, j, num;

  i=(int) floor(x);
  j=(int) floor(y);

  num=j+(i*ly);

  return num;
}

void init_frame(int lx, int ly, int *frame)
{
  int i, j, i2, j2;
  int nx, ny, n, nn, num ;
  n=0;
  
  for (i=0 ; i<lx ; i++)
    {
      for (j=0 ; j<ly ; j++)
	{
	  nn = 9*n;
	  for (i2=-1 ; i2<2 ; i2++)
	    {
	      nx = i+i2 ;
	      if (nx < 0)
		nx += lx ;
	      if (nx >= lx)
		nx -= lx;
	      nx *= ly ;
		  
	      for (j2=-1 ; j2<2 ; j2++)
		{
		  ny = j+j2 ;
		  if (ny < 0)
		    ny += ly;
		  if (ny >= ly)
		    ny -= ly ;		    		      
		    
		  num = nx+ny;
		  frame[nn++]=num;
		}/* end j2 */
	    }/* end i2 */
	  n++;
	}/* end j */
    }/* end i */
  return ;
}

void shuffleinbox2d(int lx, int ly, int nped,
		    double *x, double *y,
		    int *box, int *box0, int *neighbors, int *indx)
{
  int n, i, num, tmp0;
	
  for(num=0 ; num<(lx*ly)+1 ; num++)
    box[num]=0;

  for(n=0 ; n<nped ; n++)
    {
      num=WhichBox(ly, x[n], y[n])+1;
      neighbors[n]=num;  // neighbors[n] contains the right-neighbors box number of particle n
      box[num]+=1; // box[num] contains how many particles are left-neighbors of box number num
    }

  box0[0]=box[0];
	
  tmp0=box[1];
  box0[1]=tmp0;
  for(num=2 ; num<(lx*ly)+1 ; num++)
    {
      tmp0 += box[num] ;
      box0[num]= tmp0;
      box[num]= tmp0;
    }

  for(n=0 ; n<nped ; n++)
    {
      num=neighbors[n]-1;
      i=box0[num];
      indx[i]=n;
      box0[num]=i+1;
    }  
  return;
}


void computeDisplacement(int lx, int ly, int N, double w0, double dt, double D, double r0, double v0,
			 double beta, pcg32_random_t *rng, double *x, double *y, int *s, int *new_s,
			 int *frame, int *box, int *indx, double *dispx, double *dispy)
{
  // KEEP THREADS ALIVE ? doc openmp

#pragma omp parallel for default(shared)
  for(int i=0 ; i<N ; i++)
    {

      double xi=x[i];
      double yi=y[i];
      double si=s[i];
      int loc_m=0;

      /* FINDING ITS NEIGHBORS */
      int nn=0;
      int num=WhichBox(ly, xi, yi); 
	  
      /* Walk through neighbouring boxes */
      for (int i0=0 ; i0<9 ; i0++) 
	{
	  int num2=frame[9*num+i0];
	  int nstart=box[num2];
	  int nend=box[num2+1];
				
	  for(int ii=nstart ; ii<nend ; ii++)
	    {
	      int j=indx[ii];
	      double x1=x[j];
	      double y1=y[j];
	      double dx=x1-xi; // Vecteur orienté vers le voisin
	      double dy=y1-yi;
	      double dx1=lx+x1-xi;
	      double dy1=ly+y1-yi;
	      double dx2=-lx+x1-xi;
	      double dy2=-ly+y1-yi;

	      double absdx=fabs(dx);
	      double absdy=fabs(dy);
	      
	      if(fabs(dx1)<absdx)
		dx=dx1;
	      if(fabs(dy1)<absdy)
		dy=dy1;
	      if(fabs(dx2)<absdx)
		dx=dx2;
	      if(fabs(dy2)<absdy)
		dy=dy2;

	      double dr=sqrt(dx*dx+dy*dy);

	      if(dr<r0) /* if the neighbor j is in the interaction range */
		{ 	
		  loc_m+=s[j];
		  nn++;
		}
	    }
	}
	  
      /* Flip proba */
      double W=w0*dt*exp(-si*beta*((double)loc_m)/((double) nn));
	  
      /* New magnetization */
      double p=genrand32_real2(rng);
      if(p<W)
	si=-si;

      double etax, etay;
      BoxMuller(rng, &etax, &etay);
      
      dispx[i]=si*v0*dt + sqrt(2*D*dt)*etax;
      dispy[i]=sqrt(2*D*dt)*etay;
      new_s[i]=si;
    }

  return;
}


int updatePositions(int lx, int ly, int N, double *x, double *y, double *dispx, double *dispy,
		    int *s, int *new_s, int *N_lx, int *N_ly)
{
  int tot_mag=0; 

#pragma omp parallel for default(shared) reduction(+: tot_mag)
  for(int i=0 ; i<N ; i++)		
    {
      double xi=x[i];
      double yi=y[i];
      double si=new_s[i];

      xi += dispx[i];
      yi += dispy[i];

      /* PBC in 2D */
      if(xi<0)
	{
	  xi+=lx;
	  N_lx[i]--;
	}
      if(xi>=lx)
	{
	  xi-=lx;
	  N_lx[i]++;
	}
      if(yi<0)
	{
	  yi+=ly;
	  N_ly[i]--;
	}
      if(yi>=ly)
	{
	  yi-=ly;
	  N_ly[i]++;
	}
			
      x[i]=xi;
      y[i]=yi;
      s[i]=si;
	  
      tot_mag+=si;
    }
  return tot_mag;
}

/* RANDOM TOOLS */

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

void BoxMuller(pcg32_random_t *rng, double *x, double *y)
{
  double norm, spd2, spdx, spdy;
  
  spd2=2.;
  while ((spd2>1.)||(spd2<DBL_EPSILON))
    {
      spdx=2.*genrand32_real1(rng)-1.;
      spdy=2.*genrand32_real1(rng)-1.;

      spd2=spdx*spdx+spdy*spdy;
    }
  norm=sqrt(-2.*log(spd2)/spd2);
      
  *x = spdx*norm;
  *y = spdy*norm;

  return ;
}
