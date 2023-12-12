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
  double phi, rhol, rhog;

  /* ----------- Local observables -------- */

  int **local_m, **local_rho;
  
  /******************* TMP & INCREMENT **************/
	
  int j, i, max_mag, x_i0, y_i0;
  double clap, p, tot_mag, t1, avg_rho;
  
  /******************* Files ************************/
  	
  FILE *f_input, *f_gnuplot, *f_profiles, *f_mag, *f_plotmag;
  FILE *f_bands_rho, *f_bands_m, *f_td;
  
  f_input=fopen("f_input.dat", "r");

  f_td=fopen("f_td.dat", "w");

  f_gnuplot=fopen("plot_profiles.gp", "w");
  f_profiles=fopen("f_profiles.dat", "w");

  f_mag=fopen("f_mag.dat", "w");
  f_plotmag=fopen("plot_mag.gp", "w");

  f_bands_rho=fopen("f_bands_rho.dat", "w");
  f_bands_m=fopen("f_bands_m.dat", "w");
  
  /****************** Random initialization **********/
	
  pcg32_random_t *rng; 
  
  /* ----------------------------------------------------------------------------- */			
  /***************************************************/	
  /****************** Initialization *****************/
  /***************************************************/
  
  /****************** Parameter inputs ***************/
	
  if(fscanf(f_input, "tgap = %lg tmax = %lg rho0 = %lg lx = %d ly = %d w0 = %lg beta = %lg v = %lg D = %lg phi = %lg rhol = %lg rhog = %lg", &tgap, &tmax, &rho0, &lx, &ly, &w0, &beta, &v, &D, &phi, &rhol, &rhog)!=12){exit(1);};

  /************* Rescale parameters by detection surface **************/

  rhol=-1.35*log(1.*D)+3.8;
  
  /* D=9.*D; */
  /* phi=3*phi; */
  /* rho0=rho0/9.; */
  /* rhol=rhol/9.; */
  /* rhog=rhog/9.; */
  /* lx=3*lx; */
  /* ly=3*ly; */
  /* v=3*v; */
  
  //  N = (int) floor(rho0*lx*ly);
  
  dt=1./(4.*D+v+w0*exp(beta));

  /* int Nl = (int) floor(rhol*phi*ly); */
  /* int Ng = (int) floor(rhog*(lx-phi)*ly); */

  /* N = Nl + Ng; */

  N=(int) floor(rho0*lx*ly);
  
  printf("dt = %lg tgap = %lg tmax = %lg rhol = %lg N = %d lx = %d ly = %d w0 = %lg beta = %lg v = %lg D = %lg\n" , dt, tgap, tmax, rhol, N, lx, ly, w0, beta, v, D);
  
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
  /****************** vitesses, spins *******************/
  /******************************************************/
  
  constant();

  /* Compute Nd destabilization occurring times */

  int Nd=500;
  for(int it_Nd=0 ; it_Nd<Nd ; it_Nd++)
    {
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
  
      for(i=0 ; i<N ; i++) //for(i=0 ; i<Nl ; i++)
	{
	  x_i0=(int) floor(phi*lx*genrand32_real2(rng));
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

      /* for(i=Nl ; i<Nl+Ng ; i++) */
      /* 	{ */
      /* 	  x_i0=(int) floor(phi+(lx-phi)*genrand32_real2(rng)); */
      /* 	  x[i]=x_i0; */

      /* 	  y_i0=(int) floor(ly*genrand32_real2(rng)); */
      /* 	  y[i]=y_i0; */
      
      /* 	  p=genrand32_real2(rng); */
      /* 	  s[i]=1; */
      /* 	  if(p>.5) */
      /* 	    s[i]=-1; */

      /* 	  local_m[y_i0][x_i0]+=s[i]; */
      /* 	  local_rho[y_i0][x_i0]+=1; */
      /* 	  tot_mag+=s[i]; */
      /* 	} */
  
      t=0.;
      clap=0.;
      t1=0;
      double avg_m=0;
      double row_mag;
      double row_mag_min;
  
      while(1) /* TEMPORAL LOOP */
	{
	  /* SHUFFLE IN BOX */ 
	  tot_mag=computeLocalQuantities(lx, ly, N, x, y, s, local_m, local_rho);

	  /* COMPUTE DISPLACEMENT */
	  updatePositions(lx, ly, x, y, s, N, local_m, local_rho, w0, beta, v, D, dt, rng);

	  /* Output bands */
	  /* if(t>=clap) */
	  /*   { */
	  /*     for(i=0 ; i<lx ; i++) */
	  /* 	{ */
	  /* 	  avg_rho=0; */
	  /* 	  avg_m=0; */
	  /* 	  for(j=0 ; j<ly ; j++) */
	  /* 	    { */
	  /* 	      avg_rho+=local_rho[j][i]; */
	  /* 	      avg_m+=local_m[j][i]; */
	  /* 	    } */
	      
	  /* 	  fprintf(f_bands_rho, "%f ", ((double) avg_rho)/((double) ly)); */
	  /* 	  fprintf(f_bands_m, "%f ", ((double) avg_m)/((double) ly)); */
	  /* 	} */
	  /*     fprintf(f_bands_rho, "\n"); */
	  /*     fprintf(f_bands_m, "\n"); */
	  /*     clap+=tgap; */
	  /*   } */
      
	  /* Save density, magnetization */
	  if(t>=t1)
	    {
	      t1+=100;
	      fprintf(f_mag, "%f %f\n", t, tot_mag);
	    }

	  /* Detect destab & stop simulation right after */
	  row_mag_min=9000;
      
	  for(j=0 ; j<ly ; j++)
	    {
	      row_mag=0;
	      for (i=0 ; i<lx ; i++)
		{
		  if(fabs(local_m[j][i])>max_mag)
		    max_mag=fabs(local_m[j][i]);
		  row_mag+=local_m[j][i];
		}
	      if(row_mag<row_mag_min)
		row_mag_min=row_mag;
	    }
	  
	  if(row_mag_min<0)
	    {
	      for(i=0 ; i<lx ; i++)
		{
		  avg_rho=0;
		  avg_m=0;
		  for(j=0 ; j<ly ; j++)
		    {
		      avg_rho+=local_rho[j][i];
		      avg_m+=local_m[j][i];
		    }
	      
		  /* fprintf(f_bands_rho, "%f ", ((double) avg_rho)/((double) ly)); */
		  /* fprintf(f_bands_m, "%f ", ((double) avg_m)/((double) ly)); */
		}
	      /* fprintf(f_bands_rho, "\n"); */
	      /* fprintf(f_bands_m, "\n"); */

	      /* Save magnetization profile */
	      if(it_Nd<10)
		{
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
		}
	      
	      break;
	    }
      
	  t+=dt;
      
	} /* END WHILE */

      fprintf(f_td, "%.2f\n", t);
      fflush(f_td);
    }
	    
  /* fprintf(f_gnuplot, "set pm3d map\nset cbrange[-%f:%f]\nset palette rgb 33,13,10\nset terminal pngcairo size 800,800 enhanced font 'Verdana,10'\ndo for [i=0:%d] {\nset output sprintf('snaps/snap_%%05.0f.png', i);\n\t splot 'f_profiles.dat' index i  matrix with image notitle,\n \t set title sprintf('T = %%d', i)\n}",  rhol+1, rhol+1, (int) (tmax/tgap)-1); */

  /* fprintf(f_plotmag, "set yrange[-1.1:1.1]; set xrange[0:%d] ; set terminal png size 2000,800; set output 'total_mag.png'; plot 'f_mag.dat'", (int) tmax); */


  /* ----------------------------------------------------------------------*/	
  /******************* CLOSE FILES ************************/

  fclose(f_input);

  fclose(f_td);
  
  fclose(f_profiles);
  fclose(f_gnuplot);

  fclose(f_mag);
  fclose(f_plotmag);

  fclose(f_bands_rho);
  fclose(f_bands_m);
  
  /* if(system("gnuplot plot_mag.gp ; rm -rf snaps ; mkdir snaps ; gnuplot plot_profiles.gp ; cd snaps")!=1) {exit(1);}; */
    
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
