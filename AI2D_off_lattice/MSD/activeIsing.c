#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <float.h>
#include <string.h>
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

  int *s, *new_s;
	
  double *x, *y, *dispx, *dispy;
	
  double rho0, t, tgap, tmax, dt, r0;
  double v0, D, beta, w0;

  /* --------------- MSD ------------- */

  int *x0, *y0, *N_lx, *N_ly;
  double MSD;
  
  /* ----------------Boxes---------------- */	
  
  int *frame, *box0, *box, *indx, *neighbors;

  /******************* TMP & INCREMENT **************/
	
  int i, surf, surf1, tot_mag;
  int num, nstart, nend, ii, jj, j;
  double tmp, tmp1, clap;
  int z, itstep, min, sec, loc_mag;
  double currstep;
  clock_t t0; 
  double avg_prof, avg_rho;
    
  /******************* Files ************************/
  	
  FILE *f_coord, *f_input, *f_gp_prof, *f_profiles, *f_mag, *f_plotmag;
  FILE *f_gp_parts, *f_gp_avg_prof, *f_MSD;

  f_input=fopen("f_input.dat", "r");
  f_coord=fopen("f_coord.dat", "w");
  f_gp_parts=fopen("plot_parts.gp", "w");
  
  f_gp_prof=fopen("plot_profiles.gp", "w");
  f_profiles=fopen("f_profiles.dat", "w");

  f_MSD=fopen("f_MSD.dat", "w");
  
  f_mag=fopen("f_mag.dat", "w");
  f_plotmag=fopen("plot_mag.gp", "w");
        	
  /****************** Random initialization **********/
	
  pcg32_random_t *rng; 
  
  /* ----------------------------------------------------------------------------- */			
  /***************************************************/	
  /****************** Initialization *****************/
  /***************************************************/
  
  /****************** Parameter inputs ***************/
	
  int fscan =  fscanf(f_input, "tgap = %lg tmax = %lg rho0 = %lg lx = %d ly = %d r0 = %lg v0 = %lg D = %lg beta = %lg w0 = %lg" , &tgap, &tmax, &rho0, &lx, &ly, &r0, &v0, &D, &beta, &w0);
    
  N= (int) floor(rho0*lx*ly);
  dt = fmin((0.2*r0)*(0.2*r0)/(2.*D), (0.2*r0)/v0);
  if(dt>exp(-beta))
    printf("Decrease dt.");
    
  printf("dt = %lg tgap = %lg tmax = %lg N = %d lx = %d ly = %d r0 = %lg v0 = %lg D = %lg beta = %lg w0 = %lg\n", dt, tgap, tmax, N, lx, ly, r0, v0, D, beta, w0);

  surf=lx*ly;
  surf1=surf+1;

  char path[200];
  int pos = 0;
  pos += sprintf(&path[pos], "run_rho%.2f_", rho0);
  pos += sprintf(&path[pos], "v%.2f_", v0);
  pos += sprintf(&path[pos], "D%.2f_", D);
  pos += sprintf(&path[pos], "beta%.2f", beta);
  
  /* ----------------------------------------------------------------------------- */		
  /***************************************************/	
  /****************** Arrays allocation **************/
  /***************************************************/
	
  x=(double *)malloc(N*sizeof(double));
  EXIT_IF(x == NULL, err_malloc);
  	
  y=(double *)malloc(N*sizeof(double));
  EXIT_IF(y == NULL, err_malloc);

  s=(int *)malloc(N*sizeof(int));
  EXIT_IF(s == NULL, err_malloc);

  new_s=(int *)malloc(N*sizeof(int));
  EXIT_IF(s == NULL, err_malloc);
  
  dispx=(double *)malloc(N*sizeof(double));
  EXIT_IF(dispx == NULL, err_malloc);

  dispy=(double *)malloc(N*sizeof(double));
  EXIT_IF(dispy == NULL, err_malloc);

  /* ---------------- MSD ----------------- */

  x0=(int *)malloc(N*sizeof(int));
  EXIT_IF(x0 == NULL, err_malloc);

  y0=(int *)malloc(N*sizeof(int));
  EXIT_IF(y0 == NULL, err_malloc);

  N_lx=(int *)malloc(N*sizeof(int));
  EXIT_IF(N_lx == NULL, err_malloc);

  N_ly=(int *)malloc(N*sizeof(int));
  EXIT_IF(N_ly == NULL, err_malloc);

  /* ----------------Boxes---------------- */
  
  frame=(int *)malloc(9*surf*sizeof(int));
  EXIT_IF(frame == NULL, err_malloc);
  
  indx=(int *)malloc(N*sizeof(int));
  EXIT_IF(indx == NULL, err_malloc);
  
  box=(int *)malloc(surf1*sizeof(int));
  EXIT_IF(box == NULL, err_malloc); 	
  
  box0=(int *)malloc(surf1*sizeof(int));
  EXIT_IF(box0 == NULL, err_malloc); 
  
  neighbors=(int *)malloc(N*sizeof(int));
  EXIT_IF(neighbors == NULL, err_malloc);
  
  /****************** Seed of random number generator ***/
 	
  rng = (pcg32_random_t *)malloc(1*sizeof(pcg32_random_t));
  pcg32_srandom_r(&rng[0], time(NULL), (intptr_t)&rng[0]);

  /* ----------------------------------------------------------------------------- */	
  /******************************************************/
  /****************** Initialisation position, **********/
  /****************** vitesses, spins********************/
  /******************************************************/
  
  constant();
  t=0.;
  
  for(i=0 ; i<N ; i++)
    {
      x[i]=0.2*lx*genrand32_real2(rng);
      y[i]=ly*genrand32_real2(rng);
      x0[i]=x[i];
      y0[i]=y[i];

      N_lx[i]=0;
      N_ly[i]=0;
      
      dispx[i]=0;
      dispy[i]=0;
      s[i]=1;
      new_s[i]=1;
    }

  clap=0.;
  avg_rho=0;
  init_frame(lx, ly, frame);
  tot_mag=0;
  double t1=0;
  /* ----------------------------------------------------------------------------- */	
  /****************** Boucle dynamique ******************/	
  itstep=1;

  while(t<tmax) /* TEMPORAL LOOP */
    {
      /* SHUFFLE IN BOX */ 
      shuffleinbox2d(lx, ly, N, x, y, box, box0, neighbors, indx);

      /* COMPUTE DISPLACEMENT */
      computeDisplacement(lx, ly, N, w0, dt, D, r0, v0, beta, rng, x, y, s, new_s,
			  frame, box, indx, dispx, dispy);

      /* UPDATE POSITIONS */
      tot_mag = updatePositions(lx, ly, N, x, y, dispx, dispy, s, new_s, N_lx, N_ly);

      /* COMPUTE MSD */
      MSD=0;
      for(i=0 ; i<N ; i++)
	{
	  MSD+=(x[i]+lx*N_lx[i]-x0[i])*(x[i]+lx*N_lx[i]-x0[i])+(y[i]+ly*N_ly[i]-y0[i])*(y[i]+ly*N_ly[i]-y0[i]);
	}

      fprintf(f_MSD, "%f %f\n", t, ((double) MSD)/N);
      
      if(t>=t1)
	{
	  /* Save density and magnetization */
	  fprintf(f_mag, "%f %f\n", t, ((double) tot_mag)/N);
	  t1+=50;
	}

      
      t+=dt;		
		
    } /* END WHILE */
    
  fprintf(f_gp_prof, "set pm3d map\nset cbrange[-%d:%d]\nset palette rgb 33,13,10\nset terminal pngcairo size 800,800 enhanced font 'Verdana,10'\ndo for [i=0:%d] {\nset output sprintf('%s/snaps_prof/snap_%%05.0f.png', i);\n\t splot 'f_profiles.dat' index i  matrix with image notitle,\n \t set title sprintf('T = %%d', i)\n}", 10, 10, (int) (tmax/tgap)-1, path);

  /* fprintf(f_gp_parts, "set size ratio 1\nset xrange [0:%d]\nset yrange [0:%d]\nset palette model RGB defined (-1 'red', 1 'green' )\nunset key\nunset colorbox\nset terminal pngcairo size 800,800 enhanced font 'Verdana,10'\na=0\nb=%d\nincr=1\ni=0\ndo for [ii=a:b:incr] {\nset output sprintf('%s/snaps_parts/snaps_%%05.0f.png', i)\ni=i+1\nplot 'f_coord.dat' u 2:3:4 every :incr::ii::ii w p pt 7 palette\n}", lx, ly, (int) (tmax/tgap)-1, path); */

    /* fprintf(f_gp_avg_prof, "set size ratio 1\nset xrange [0:%d]\nset yrange [-5:5]\nunset key\nset terminal pngcairo size 800,800 enhanced font 'Verdana,10'\ndo for [i=0:%d] {\nset output sprintf('%s/snaps_avg/snaps_%%05.0f.png', i)\nplot 'f_avg_prof.dat' every :::i::i w l\n}", lx, (int) (tmax/tgap)-1, path); */

  fprintf(f_plotmag, "set yrange[-1.1:1.1]; set xrange[0:%d] ; set terminal png size 800,800; set output '%s/total_mag.png'; plot 'f_mag.dat'", (int) tmax, path);

  /* -----------------------------------------------------------------------------*/	
  /******************* CLOSE FILES ************************/

  fclose(f_input);
  fclose(f_coord);
  fclose(f_gp_parts);
  
  fclose(f_profiles);
  fclose(f_gp_prof);

  fclose(f_MSD);

  fclose(f_mag);
  fclose(f_plotmag);

  /* -----------------------------------------------------------------------------*/
  /****************** FILMS  *************************/
  
  char cmd[3000];

  strcpy(cmd, "rm -rf ");
  strcat(cmd, path);
  
  strcat(cmd, " ; mkdir ");
  strcat(cmd, path);

  strcat(cmd, " ; mkdir ");
  strcat(cmd, path);
  strcat(cmd, "/snaps_prof");

  /* strcat(cmd, " ; mkdir "); */
  /* strcat(cmd, path); */
  /* strcat(cmd, "/snaps_parts"); */

  strcat(cmd, " ; mkdir ");
  strcat(cmd, path);
  strcat(cmd, "/snaps_avg");

  /* strcat(cmd, " ; gnuplot plot_mag.gp ; gnuplot plot_profiles.gp ; gnuplot plot_parts.gp ; gnuplot plot_avg_profiles.gp"); */

  strcat(cmd, " ; gnuplot plot_mag.gp ; gnuplot plot_profiles.gp ; gnuplot plot_avg_profiles.gp");

  fscan=system(cmd);

  /* -----------------------------------------------------------------------------*/
  /****************** FREE MEMORY *************************/
  
  free(x);
  free(y);
  
  free(rng);
  
  free(frame);
  free(indx);
  free(box0);
  free(box);
  free(neighbors);
  
  free(s);
  free(new_s);
  
  return EXIT_SUCCESS;
}
