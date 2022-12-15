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

int main(void)
{
  /**************************************************/
  /******************* Declaration ******************/
  /**************************************************/

  /********* Global parameters and dynamics *********/
  
  int N, lx, ly;

  int *s, *x, *y;
	
  double rho0, dt, t, tgap, tmax;
  double v, beta, w0, D;

  /* ----------- Local observables -------- */

  int **local_m, **local_rho;
  
  /******************* TMP & INCREMENT **************/
	
  int j, i, max_mag, x_i0, y_i0;
  double clap, p, tot_mag, t1, avg_rho;
  
  /******************* Files ************************/
  	
  FILE *f_input, *f_profiles;
  
  f_input=fopen("f_input.dat", "r");
  f_profiles=fopen("f_profiles.dat", "w");
  
  /****************** Random initialization **********/
	
  pcg32_random_t *rng; 
  
  /* ----------------------------------------------------------------------------- */			
  /***************************************************/	
  /****************** Initialization *****************/
  /***************************************************/
  
  /****************** Parameter inputs ***************/
	
  if(fscanf(f_input, "tgap = %lg tmax = %lg rho0 = %lg lx = %d ly = %d w0 = %lg beta = %lg v = %lg D = %lg", &tgap, &tmax, &rho0, &lx, &ly, &w0, &beta, &v, &D)!=9){exit(1);};

  
  N = (int) floor(rho0*lx*ly);

  dt=1./(4.*D+v+w0*exp(beta));
  
  /* ----------------------------------------------------------------------------- */		
  /***************************************************/	
  /****************** Arrays allocation **************/
  /***************************************************/
	
  x=(int *)malloc((N)*sizeof(int));
  EXIT_IF(x == NULL, err_malloc);

  y=(int *)malloc((N)*sizeof(int));
  EXIT_IF(x == NULL, err_malloc);
  
  s=(int *)malloc((N)*sizeof(int));
  EXIT_IF(s == NULL, err_malloc);

  /* ------------- Observables -------------- */

  local_m=imatrix((long) ly, (long) lx);
  local_rho=imatrix((long) ly, (long) lx);
  
  /****************** Seed of random number generator ***/
 	
  rng = (pcg32_random_t *)malloc(1*sizeof(pcg32_random_t));
  pcg32_srandom_r(&rng[0], time(NULL), (intptr_t)&rng[0]);

  /* ---------------------------------------------------------------------- */	
  /******************************************************/
  /****************** Initialisation position, **********/
  /****************** vitesses, spins *******************/
  /******************************************************/
  
  constant();

  for(int it=0 ; it<1; it++)
    {
      tot_mag=0;
  
      for(j=0 ; j<ly ; j++)
	{
	  for(i=0 ; i<lx ; i++)
	    {
	      local_m[j][i]=0;
	      local_rho[j][i]=0;
	    }
	}
  
      for(i=0 ; i<N ; i++)
	{
	  x_i0=(int) floor(lx*genrand32_real2(rng));
	  x[i]=x_i0;

	  y_i0=(int) floor(ly*genrand32_real2(rng));
	  y[i]=y_i0;
      
	  p=genrand32_real2(rng);
	  s[i]=1;

	  local_m[y_i0][x_i0]+=s[i];
	  local_rho[y_i0][x_i0]+=1;
	  tot_mag+=s[i];
	}
  
      t=0.;
      clap=100;

      while(t<tmax) /* TEMPORAL LOOP */
	{
      
	  /* SHUFFLE IN BOX */ 
	  double frac_rev=computeLocalQuantities(lx, ly, N, x, y, s, local_m, local_rho); 
      
	  /* COMPUTE DISPLACEMENT */
	  updatePositions(lx, ly, x, y, s, N, local_m, local_rho, w0, beta, v, D, dt, rng);

	  t+=dt;

	  if(t>=clap)
	    {
	      /* Save magnetization profile */
	      for(j=0 ; j<ly ; j++)
		{
		  for (i=0 ; i<lx ; i++)
		    {
		      fprintf(f_profiles, "%f ", ((double)local_m[j][i]));
		    }
		  fprintf(f_profiles, "\n");
		}
	      fprintf(f_profiles, "\n\n");
	      fflush(f_profiles);
	      clap += tgap;
	    }      
	} /* END WHILE */
    }
  /* ----------------------------------------------------------------------*/	
  /******************* CLOSE FILES ************************/

  fclose(f_input);

  
  fclose(f_profiles);
    
  /* -----------------------------------------------------------------------------*/
  /****************** FREE MEMORY *************************/
	
  free(x);
  free(y);
  
  free(rng);
  
  free_imatrix(local_m);
  free_imatrix(local_rho);
  
  free(s);
  
  return EXIT_SUCCESS;
}
