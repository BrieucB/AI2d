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
  double eps, beta, w0, D;

  /* ----------- Local observables -------- */

  int **local_m, **local_rho;
  
  /******************* TMP & INCREMENT **************/
	
  int j, i, max_mag, x_i0, y_i0;
  double clap, p, tot_mag, t1, avg_rho;
  
  /******************* Files ************************/
  	
  FILE *f_input, *f_gnuplot, *f_mag, *f_plotmag;
  FILE *f_profiles_rho, *f_profiles_m;
  
  f_input=fopen("f_input.dat", "r");

  f_gnuplot=fopen("plot_profiles.gp", "w");

  f_mag=fopen("f_mag.dat", "w");
  f_plotmag=fopen("plot_mag.gp", "w");

  f_profiles_rho=fopen("f_profiles_rho.dat", "w");
  f_profiles_m=fopen("f_profiles_m.dat", "w");
  
  /****************** Random initialization **********/
	
  pcg32_random_t *rng; 
  
  /* ----------------------------------------------------------------------------- */			
  /***************************************************/	
  /****************** Initialization *****************/
  /***************************************************/
  
  /****************** Parameter inputs ***************/
	
  if(fscanf(f_input, "tgap = %lg tmax = %lg rho0 = %lg lx = %d ly = %d w0 = %lg beta = %lg eps = %lg D = %lg", &tgap, &tmax, &rho0, &lx, &ly, &w0, &beta, &eps, &D)!=9){exit(1);};

  D=9.*D;
  rho0=rho0/9.;
  lx=3*lx;
  ly=3*ly;
    
  N= (int) floor(rho0*lx*ly);
  dt=1./(4.*D+w0*exp(beta));
  
  printf("dt = %lg tgap = %lg tmax = %lg rho0 = %lg N = %d lx = %d ly = %d w0 = %lg beta = %lg eps = %lg D = %lg\n" , dt, tgap, tmax, rho0, N, lx, ly, w0, beta, eps, D);

  
  /* ----------------------------------------------------------------------------- */		
  /***************************************************/	
  /****************** Arrays allocation **************/
  /***************************************************/
	
  x=(int *)malloc(N*sizeof(int));
  EXIT_IF(x == NULL, err_malloc);

  y=(int *)malloc(N*sizeof(int));
  EXIT_IF(x == NULL, err_malloc);
  
  s=(int *)malloc(N*sizeof(int));
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
  /****************** vitesses, spins********************/
  /******************************************************/
  
  constant();
  max_mag=0;
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
      x_i0=(int) floor(0.8*lx*genrand32_real2(rng));
      x[i]=x_i0;

      y_i0=(int) floor(ly*genrand32_real2(rng));
      y[i]=y_i0;
      
      p=genrand32_real2(rng);
      s[i]=1;
      if(p>1.5)
	s[i]=-1;

      local_m[y_i0][x_i0]+=s[i];
      local_rho[y_i0][x_i0]+=1;
      tot_mag+=s[i];
    }

  t=0.;
  clap=0.;
  t1=0;

  double t_destab=99999999999;
  int not_destab=1;
  double t_reset=0;
  
  /* ---------------------------------------------------------------------- */	
  /****************** Boucle dynamique ******************/	
  
  while(t<tmax) /* TEMPORAL LOOP */
    {
      /* COMPUTE FIELDS */ 
      tot_mag=computeLocalQuantities(lx, ly, N, x, y, s, local_m, local_rho);

      /* COMPUTE DISPLACEMENT */
      updatePositions(lx, ly, x, y, s, N, local_m, local_rho, w0, beta, eps, D, dt, rng);

      if(t>=clap)
	{
	  if(t>=t_reset)
	    {
	      fseek(f_profiles_m, 0, SEEK_SET);
	      fseek(f_profiles_rho, 0, SEEK_SET);
	      t_reset+=20000;
	    }

	  /* Output magnetization in each box */
	  for(j=0 ; j<ly ; j++)
	    {

	      for(i=0 ; i<lx ; i++)
		{

		  if(i==0)
		    {
		      fprintf(f_profiles_m, "%.2f %d ", t, local_m[j][i]);
		      fprintf(f_profiles_rho, "%.2f %d ", t, local_rho[j][i]);
		    }

		  else
		    {
		      fprintf(f_profiles_m, "%d ", local_m[j][i]);
		      fprintf(f_profiles_rho, "%d ", local_rho[j][i]);
		    }

		}

	      fprintf(f_profiles_m, "\n");
	      fprintf(f_profiles_rho, "\n");	
		 
	    }

	  fprintf(f_profiles_m, "\n\n");
	  fprintf(f_profiles_rho, "\n\n");

	  if((tot_mag<0.)&&(not_destab==1))
	    {
	      t_destab=t+10000;
	      not_destab=0;	    
	    }
	      
	  clap+=tgap;
	}

      if(t>=t1)
	{
	  t1+=100;
	  fprintf(f_mag, "%f %f\n", t, tot_mag);
	  fflush(f_mag);
	}
      
      /* /\* Save density, magnetization *\/ */
      /* if(t>=t1) */
      /* 	{ */
      /* 	  t1+=10; */
      /* 	  fprintf(f_mag, "%f %f\n", t, tot_mag); */
      /* 	  fflush(f_mag); */
      /* 	} */
      
      /* if(t>=clap) */
      /* 	{ */
      /* 	  for(j=0 ; j<ly ; j++) */
      /* 	    { */
      /* 	      for (i=0 ; i<lx ; i++) */
      /* 		{ */
      /* 		  if(i==0) */
      /* 		    { */
      /* 		      fprintf(f_profiles, "%.2f %d ", t, local_m[j][i]); */
      /* 		    } */
		  
      /* 		  else */
      /* 		    { */
      /* 		      fprintf(f_profiles, "%d ", local_m[j][i]); */
      /* 		    } */
      /* 		} */

      /* 	      fprintf(f_profiles, "\n"); */
      /* 	    } */

      /* 	  fprintf(f_profiles, "\n\n"); */
      /* 	  clap+=tgap; */
      /* 	} */
      
      t+=dt;
      
    } /* END WHILE */
	    
  fprintf(f_gnuplot, "set pm3d map\nset cbrange[-%d:%d]\nset palette rgb 33,13,10\nset terminal pngcairo size 800,800 enhanced font 'Verdana,10'\ndo for [i=0:%d] {\nset output sprintf('snaps/snap_%%05.0f.png', i);\n\t splot 'f_profiles.dat' index i  matrix with image notitle,\n \t set title sprintf('T = %%d', i)\n}", max_mag+1, max_mag+1, (int) (tmax/tgap)-1);

  fprintf(f_plotmag, "set yrange[-1.1:1.1]; set xrange[0:%d] ; set terminal png size 2000,800; set output 'total_mag.png'; plot 'f_mag.dat'", (int) tmax);


  /* ----------------------------------------------------------------------*/	
  /******************* CLOSE FILES ************************/

  fclose(f_input);

  fclose(f_gnuplot);

  fclose(f_mag);
  fclose(f_plotmag);

  fclose(f_profiles_rho);
  fclose(f_profiles_m);
  
  /* if(system("gnuplot plot_mag.gp ; rm -rf snaps ; mkdir snaps ; gnuplot plot_profiles.gp ; cd snaps")!=1) {exit(1);}; */
    
  /* -----------------------------------------------------------------------------*/
  /****************** FREE MEMORY *************************/
	
  free(x);
  
  free(rng);
  
  free_imatrix(local_m);
  free_imatrix(local_rho);
  
  free(s);
  
  return EXIT_SUCCESS;
}
