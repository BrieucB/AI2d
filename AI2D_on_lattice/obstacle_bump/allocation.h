#ifndef  ALLOCATION_H
#define  ALLOCATION_H

/* fonctions */
double **dmatrix(long , long);
int **imatrix(long , long);
double ***dtensor3(long, long, long);
int ***itensor3(long, long, long);
double ****dtensor4(long, long, long, long);
int ****itensor4(long, long, long, long);
void free_dmatrix(double**);
void free_imatrix(int**);
void free_dtensor3(double***);
void free_itensor3(int***);
void free_dtensor4(double****);
void free_itensor4(int****);

void re_imatrix(int **,long, long ) ;
void re_dmatrix(double **,long , long ) ;

#ifdef FFTW_IN_USE
fftw_complex **fftwmatrix(long , long );
fftw_complex ***fftwtensor3(long , long , long );
void free_fftwmatrix(fftw_complex **);
void free_fftwtensor3(fftw_complex ***);
#endif /* FFTW_IN_USE */

#endif /* ALLOCATION _H*/
