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
  double h0;
  int sigma0, Nit; 
  
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
	
  if(fscanf(f_input, "tgap = %lg tmax = %lg rho0 = %lg lx = %d ly = %d w0 = %lg beta = %lg v = %lg D = %lg h0 = %lg r0 = %d Nit = %d", &tgap, &tmax, &rho0, &lx, &ly, &w0, &beta, &v, &D, &h0, &sigma0, &Nit)!=12){exit(1);};

  
  N = (int) floor(rho0*lx*ly);

  h0*=rho0;
  int Nfluct = (int) floor(M_PI*h0*sigma0*sigma0);
  
  dt=1./(4.*D+v+w0*exp(beta));

  
  printf("dt = %lg tgap = %lg tmax = %lg N = %d lx = %d ly = %d w0 = %lg beta = %lg v = %lg D = %lg h0 = %lg r0 = %d\n" , dt, tgap, tmax, N, lx, ly, w0, beta, v, D, h0, sigma0);
  
  /* ----------------------------------------------------------------------------- */		
  /***************************************************/	
  /****************** Arrays allocation **************/
  /***************************************************/
	
  x=(int *)malloc((N+Nfluct)*sizeof(int));
  EXIT_IF(x == NULL, err_malloc);

  y=(int *)malloc((N+Nfluct)*sizeof(int));
  EXIT_IF(x == NULL, err_malloc);
  
  s=(int *)malloc((N+Nfluct)*sizeof(int));
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

  for(int it=0 ; it<Nit; it++)
    {
      N = (int) floor(rho0*lx*ly);
      int Nfluct = (int) floor(M_PI*h0*sigma0*sigma0);

      max_mag=0;
      tot_mag=0;
      double frac_rev=0.;
      double frac_rev0=0;
  
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
      clap=0.;
      t1=0;
      double tblob=10;
      double avg_m=0;
      int lx2 = (int) floor(((double)lx)/2);
      int ly2 = (int) floor(((double)ly)/2);
      int destab=0;
      
      while( (t<tblob) || (destab==0) ) /* TEMPORAL LOOP */
	{
	  /* ADD A FLUCTUATION */
	  if((t>=tblob)&&(t<tblob+dt))
	    {

	      /* By reversing the particles in the fluctuation zone */
	      for(i=0 ; i<N ; i++)
		{
		  x_i0=x[i];
		  y_i0=y[i];

		  if((lx2-x_i0)*(lx2-x_i0) + (ly2-y_i0)*(ly2-y_i0) <= sigma0*sigma0)
		    {
		      s[i]=-1;
		    }
		}
	      
	      /* By adding particles */
	      N+=Nfluct;
	      for(i=N-Nfluct ; i<N ; i++)
		{
		  x_i0=(int) floor(lx2-sigma0 + 2.*sigma0*genrand32_real2(rng));
		  y_i0=(int) floor(ly2-sigma0 + 2.*sigma0*genrand32_real2(rng));
		  
		  while((lx2-x_i0)*(lx2-x_i0) + (ly2-y_i0)*(ly2-y_i0) > sigma0*sigma0)
		    {
		      x_i0=(int) floor(lx2-sigma0 + 2.*sigma0*genrand32_real2(rng));
		      y_i0=(int) floor(ly2-sigma0 + 2.*sigma0*genrand32_real2(rng));		  
		    }

		  x[i]=x_i0;
		  y[i]=y_i0;
		  s[i]=-1;
	      
		  local_m[y_i0][x_i0]+=s[i];
		  local_rho[y_i0][x_i0]+=1;
		  tot_mag+=s[i];
		}
	    }


	  
	  /* SHUFFLE IN BOX */ 
	  frac_rev=computeLocalQuantities(lx, ly, N, x, y, s, local_m, local_rho); 
      
	  /* COMPUTE DISPLACEMENT */
	  updatePositions(lx, ly, x, y, s, N, local_m, local_rho, w0, beta, v, D, dt, rng);

	  /* /\* Save magnetization *\/ */
	  /* if(t>=t1) */
	  /*   { */
	  /*     t1+=0.1; */
	  /*     /\* fprintf(f_mag, "%f %f\n", t, tot_mag); *\/ */
	  /*     fprintf(f_mag, "%f %f\n", t, frac_rev); */
	  /*     fflush(f_mag); */
	  /*   } */

	  /* if(t>=clap) */
	  /*   { */
	  /*     /\* Save magnetization profile *\/ */
	  /*     for(j=0 ; j<ly ; j++) */
	  /* 	{ */
	  /* 	  for (i=0 ; i<lx ; i++) */
	  /* 	    { */
	  /* 	      fprintf(f_profiles, "%f ", ((double)local_m[j][i])); */
	  /* 	    } */
	  /* 	  fprintf(f_profiles, "\n"); */
	  /* 	} */
	  /*     fprintf(f_profiles, "\n\n"); */
	  /*     fflush(f_profiles); */

	  /*     fprintf(f_mag, "%f %f\n", t, frac_rev); */
	      
	  /*     clap+=tgap; */
	  /*   } */

	  int mag_ahead=0;
	  for(int j=ly2-2 ; j<ly2+3 ; j++ )
	    {
	      for(int i=lx2-50-2 ; i<lx2-50+3 ; i++ )
		{
		  mag_ahead+=local_m[j][i];
		}
	    }
	  if(mag_ahead<-25*rho0)
	    destab=1;

	  if(t>=200)
	    destab=1;
	  
	  t+=dt;
	  //	  printf("%f %f %f %f\n", t, frac_rev0, frac_rev, 0.1*100*200/(lx*ly));
	} /* END WHILE */

      if(t>=200)      
	fprintf(f_mag, "0 ");
      else
	fprintf(f_mag, "%f ", t);
      
      fflush(f_mag);
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
