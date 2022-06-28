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

  /* ----------------Boxes---------------- */	
  
  int *frame, *box0, *box, *indx, *neighbors;

  /******************* TMP & INCREMENT **************/
	
  int i, surf, surf1, tot_mag, N_run, k_run;
  int num, nstart, nend, ii, jj, j;
  double clap;
  int loc_mag;
  int not_destab, row_mag, row_mag_min;
  
  /******************* Files ************************/
  	
  FILE *f_coord, *f_input, *f_gp_prof, *f_profiles_rho, *f_profiles_m, *f_mag, *f_plotmag;
  FILE *f_gp_parts, *f_gp_avg_prof, *f_avg_prof;
  //FILE *f_init;
  
  f_input=fopen("f_input.dat", "r");
  f_coord=fopen("f_coord.dat", "w");
  f_gp_parts=fopen("plot_parts.gp", "w");
  
  f_gp_prof=fopen("plot_profiles.gp", "w");

  f_gp_avg_prof=fopen("plot_avg_profiles.gp", "w");
  f_avg_prof=fopen("f_avg_prof.dat", "w");
  
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
  dt = fmin((0.1*r0)*(0.1*r0)/(2.*D), (0.1*r0)/v0);
  if(dt>exp(-beta))
    printf("Decrease dt.");
    
  printf("dt = %lg tgap = %lg tmax = %lg N = %d lx = %d ly = %d r0 = %lg v0 = %lg D = %lg beta = %lg w0 = %lg\n", dt, tgap, tmax, N, lx, ly, r0, v0, D, beta, w0);

  /* char init_name[200]; */
  /* int pos_0 = 0; */
  /* pos_0 += sprintf(&init_name[pos_0], "pos_init_D%lg.dat", D); */
  /* f_init=fopen(init_name, "r"); */
  
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
  N_run=10;

  for(k_run=0 ; k_run<N_run ; k_run++)
    { 
     t=0.;

     char f_name_m[200];
     int pos_m = 0;
     pos_m += sprintf(&f_name_m[pos_m], "f_profiles_m_run%d.dat", k_run);

     f_profiles_m=fopen(f_name_m, "w");

     char f_name_rho[200];
     int pos_rho = 0;
     pos_rho += sprintf(&f_name_rho[pos_rho], "f_profiles_rho_run%d.dat", k_run);

     f_profiles_rho=fopen(f_name_rho, "w");

     
     
      for(i=0 ; i<N ; i++)
	{
	  /* fscan =  fscanf(f_init, "%.3f %.3f %d\n", &x[i], &y[i], &s[i]);  */
	  x[i]=0.2*lx*genrand32_real2(rng);
	  y[i]=ly*genrand32_real2(rng);
	  dispx[i]=0;
	  dispy[i]=0;
	  s[i]=1;
	  new_s[i]=s[i];
	}

      clap=0.;
      double t_reset=0;
      init_frame(lx, ly, frame);
      tot_mag=0;
      not_destab=1;
      double t_destab=999999;
      
      /* ----------------------------------------------------------------------------- */	
      /****************** Boucle dynamique ******************/	

      while(t<t_destab) /* TEMPORAL LOOP */
	{
	  /* SHUFFLE IN BOX */ 
	  shuffleinbox2d(lx, ly, N, x, y, box, box0, neighbors, indx);

	  /* COMPUTE DISPLACEMENT */
	  computeDisplacement(lx, ly, N, w0, dt, D, r0, v0, beta, rng, x, y, s, new_s,
			      frame, box, indx, dispx, dispy);

	  /* UPDATE POSITIONS */
	  tot_mag = updatePositions(lx, ly, N, x, y, dispx, dispy, s, new_s);

	  if(t>=clap)
	    {
	      if(t>=t_reset)
		{
		  fseek(f_profiles_m, 0, SEEK_SET);
		  fseek(f_profiles_rho, 0, SEEK_SET);
		  t_reset+=200;
		}
	      
	      /* Output magnetization in each box */
	      row_mag_min=0;
	      for(j=0 ; j<ly ; j++)
		{
		  row_mag=0;
		  for(i=0 ; i<lx ; i++)
		    {
		      num=i*ly+j;
		      nstart=box[num];
		      nend=box[num+1];
		      loc_mag=0;
		  
		      for(ii=nstart ; ii<nend ; ii++)
			{
			  jj=indx[ii];
			  loc_mag+=s[jj];
			}


		      if(i==0)
			{
			  fprintf(f_profiles_m, "%.2f %d ", t, loc_mag);
			  fprintf(f_profiles_rho, "%.2f %d ", t, nend-nstart);
			}

		      else
			{
			  fprintf(f_profiles_m, "%d ", loc_mag);
			  fprintf(f_profiles_rho, "%d ", nend-nstart);
			}

		      row_mag+=loc_mag;
		    }

		  if(row_mag<row_mag_min)
		    row_mag_min=row_mag;

		  fprintf(f_profiles_m, "\n");
		  fprintf(f_profiles_rho, "\n");	
		 
		}

	      fprintf(f_profiles_m, "\n\n");
	      fprintf(f_profiles_rho, "\n\n");
	      
	      if((row_mag_min<0)&&(not_destab==1))
		{
		  t_destab=t+100;
		  not_destab=0;
		  
		}

	      /* Output particles positions */
	      /* f_coord=freopen("f_coord.dat", "w", f_coord); */
	      /* for(i=0 ; i<N ; i++) */
	      /*   fprintf(f_coord, "%f %f %f %d\n", t, x[i], y[i], s[i]); */

	      /* /\* Output magnetization in each box *\/ */
	      /* for(i=0 ; i<lx ; i++) */
	      /*   { */

	      /*     avg_prof=0.; */
	      /*     avg_rho=0.; */
	      /*     for(j=0 ; j<ly ; j++) */
	      /* 	{ */
	      /* 	  num=i*ly+j; */
	      /* 	  nstart=box[num]; */
	      /* 	  nend=box[num+1]; */
	      /* 	  loc_mag=0; */
		  
	      /* 	  for(ii=nstart ; ii<nend ; ii++) */
	      /* 	    { */
	      /* 	      jj=indx[ii]; */
	      /* 	      loc_mag+=s[jj]; */
	      /* 	    } */
		  
	      /* 	  /\* Output profiles averaged along the y-axis *\/ */
	      /* 	  avg_prof+=((double)loc_mag); */
	      /* 	  avg_rho+=(double) (nend-nstart); */
	      /* 	} */
	      /*     fprintf(f_avg_prof, "%d %.2f %.2f\n", i, avg_prof/ly,avg_rho/ly); */
	      /*   } */
	      /* fprintf(f_avg_prof, "\n"); */
	  
	      /* Save total magnetization */
	      fprintf(f_mag, "%f %f\n", t, ((double) tot_mag)/N);

	      clap+=tgap;
	    }
	  
	  t+=dt;		
		
	} /* END WHILE */

    }
      
  /* fprintf(f_gp_prof, "set pm3d map\nset cbrange[-%d:%d]\nset palette rgb 33,13,10\nset terminal pngcairo size 800,800 enhanced font 'Verdana,10'\ndo for [i=0:%d] {\nset output sprintf('%s/snaps_prof/snap_%%05.0f.png', i);\n\t splot 'f_profiles.dat' index i  matrix with image notitle,\n \t set title sprintf('T = %%d', i)\n}", 10, 10, (int) (tmax/tgap)-1, path); */

  /* fprintf(f_gp_parts, "set size ratio 1\nset xrange [0:%d]\nset yrange [0:%d]\nset palette model RGB defined (-1 'red', 1 'green' )\nunset key\nunset colorbox\nset terminal pngcairo size 800,800 enhanced font 'Verdana,10'\na=0\nb=%d\nincr=1\ni=0\ndo for [ii=a:b:incr] {\nset output sprintf('%s/snaps_parts/snaps_%%05.0f.png', i)\ni=i+1\nplot 'f_coord.dat' u 2:3:4 every :incr::ii::ii w p pt 7 palette\n}", lx, ly, (int) (tmax/tgap)-1, path); */

    /* fprintf(f_gp_avg_prof, "set size ratio 1\nset xrange [0:%d]\nset yrange [-5:5]\nunset key\nset terminal pngcairo size 800,800 enhanced font 'Verdana,10'\ndo for [i=0:%d] {\nset output sprintf('%s/snaps_avg/snaps_%%05.0f.png', i)\nplot 'f_avg_prof.dat' every :::i::i w l\n}", lx, (int) (tmax/tgap)-1, path); */

  /* fprintf(f_plotmag, "set yrange[-1.1:1.1]; set xrange[0:%d] ; set terminal png size 800,800; set output '%s/total_mag.png'; plot 'f_mag.dat'", (int) tmax, path); */

  /* -----------------------------------------------------------------------------*/	
  /******************* CLOSE FILES ************************/

  fclose(f_input);
  fclose(f_coord);
  fclose(f_gp_parts);
  
  fclose(f_profiles_m);
  fclose(f_profiles_rho);
  fclose(f_gp_prof);

  fclose(f_avg_prof);
  fclose(f_gp_avg_prof);

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

  /* fscan=system(cmd); */

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
