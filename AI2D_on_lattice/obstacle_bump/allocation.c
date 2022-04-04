#include <stdlib.h>

#ifdef FFTW_IN_USE
#include <complex.h>
#include <fftw3.h>
#endif /* FFTW_IN_USE */

#include "exit_if.h"
#include "erreurs.h"
#include "allocation.h"

#define FREE_ARG char*

double **
dmatrix(long nrow, long ncol)
/* allocate a double matrix with subscript range m[0..nrow-1][0..ncol-1] */
{
  long i;
  double **m;
  
  /* allocate pointers to rows */
  m=(double **) malloc((size_t)(nrow*sizeof(double*)));
  EXIT_IF(m == NULL,err_malloc);  
  /* allocate rows and set pointers to them */
  m[0]=(double *) malloc((size_t)((nrow*ncol)*sizeof(double)));
  EXIT_IF(m[0] == NULL,err_malloc);  
  
  for(i=1;i<nrow;i++)
    m[i]=m[i-1]+ncol;
  
  /* return pointer to array of pointers to rows */
  return m;
}

void
re_dmatrix(double **m,long nrow, long ncol)
/* allocate an int matrix with subscript range m[0..nrow-1][0..ncol-1] */
{
  long i;
  
  /* allocate pointers to rows */
  m=(double **) realloc(m,(size_t)(nrow*sizeof(double *)));
  EXIT_IF(m == NULL, err_realloc);  
  
  /* allocate rows and set pointers to them */
  m[0]=(double *) realloc(m[0],(size_t)((nrow*ncol)*sizeof(double)));
  EXIT_IF(m[0] == NULL, err_realloc);  
  
  for(i=1;i<nrow;i++)
    m[i]=m[i-1]+ncol;
  
  return ;
}

void
free_dmatrix(double **m)
     /* free a double matrix allocated by dmatrix() */
{
  free((FREE_ARG) (m[0]));
  free((FREE_ARG) (m));
}

int **
imatrix(long nrow, long ncol)
/* allocate an int matrix with subscript range m[0..nrow-1][0..ncol-1] */
{
  long i;
  int **m;
  
  /* allocate pointers to rows */
  m=(int **) malloc((size_t)(nrow*sizeof(int*)));
  EXIT_IF(m == NULL,err_malloc);  
  
  /* allocate rows and set pointers to them */
  m[0]=(int *) malloc((size_t)((nrow*ncol)*sizeof(int)));
  EXIT_IF(m[0] == NULL,err_malloc);  
  
  for(i=1;i<nrow;i++)
    m[i]=m[i-1]+ncol;
  
  /* return pointer to array of pointers to rows */
  return m;
}

void
re_imatrix(int **m,long nrow, long ncol)
/* allocate an int matrix with subscript range m[0..nrow-1][0..ncol-1] */
{
  long i;
  
  /* allocate pointers to rows */
  m=(int **) realloc(m,(size_t)(nrow*sizeof(int*)));
  EXIT_IF(m == NULL,err_malloc);  
  
  /* allocate rows and set pointers to them */
  m[0]=(int *) realloc(m[0],(size_t)((nrow*ncol)*sizeof(int)));
  EXIT_IF(m[0] == NULL,err_malloc);  
  
  for(i=1;i<nrow;i++)
    m[i]=m[i-1]+ncol;
  
  /* return pointer to array of pointers to rows */
  return ;
}

void
free_imatrix(int **m)
     /* free a double matrix allocated by dmatrix() */
{
  free((FREE_ARG) (m[0]));
  free((FREE_ARG) (m));
}

double ***
dtensor3(long index1, long index2, long index3)
     /* allocate a int 3tensor */
{
  long i,j;
  double ***t;
  
  /* allocate pointers to pointers to rows */
  t=(double ***) malloc((size_t)(index1*sizeof(double**)));
  EXIT_IF(t == NULL,err_malloc);  
  
  /* allocate pointers to rows and set pointers to them */
  t[0]=(double **) malloc((size_t)((index1*index2)*sizeof(double*)));
  EXIT_IF(t[0] == NULL,err_malloc);
  
  /* allocate rows and set pointers to them */
  t[0][0]=(double *) malloc((size_t)((index1*index2*index3)*sizeof(double)));
  EXIT_IF(t[0][0] == NULL,err_malloc);  
  
  for(j=1;j<index2;j++)
    t[0][j]=t[0][j-1]+index3;
  for(i=1;i<index1;i++)
    {
      t[i]=t[i-1]+index2;
      t[i][0]=t[i-1][0]+index2*index3;
      for(j=1;j<index2;j++)
	t[i][j]=t[i][j-1]+index3;
    }

  /* return pointer to array of pointers to rows */
  return t;
}

void
free_dtensor3(double ***m)
{
  free((FREE_ARG) (m[0][0]));
  free((FREE_ARG) (m[0]));
  free((FREE_ARG) (m));
}

int ***
itensor3(long index1, long index2, long index3)
     /* allocate a int 3tensor */
{
  long i,j;
  int ***t;
  
  /* allocate pointers to pointers to rows */
  t=(int ***) malloc((size_t)(index1*sizeof(int**)));
  EXIT_IF(t == NULL,err_malloc);  
  
  /* allocate pointers to rows and set pointers to them */
  t[0]=(int **) malloc((size_t)((index1*index2)*sizeof(int*)));
  EXIT_IF(t[0] == NULL,err_malloc);  
  
  /* allocate rows and set pointers to them */
  t[0][0]=(int *) malloc((size_t)((index1*index2*index3)*sizeof(int)));
  EXIT_IF(t[0][0] == NULL,err_malloc);  
  
  for(j=1;j<index2;j++)
    t[0][j]=t[0][j-1]+index3;
  for(i=1;i<index1;i++)
    {
      t[i]=t[i-1]+index2;
      t[i][0]=t[i-1][0]+index2*index3;
      for(j=1;j<index2;j++)
	t[i][j]=t[i][j-1]+index3;
    }

  /* return pointer to array of pointers to rows */
  return t;
}

void
free_itensor3(int ***m)
{
  free((FREE_ARG) (m[0][0]));
  free((FREE_ARG) (m[0]));
  free((FREE_ARG) (m));
}

double ****
dtensor4(long index1, long index2, long index3, long index4)
     /* allocate a double 4d_tensor */
{
  long i,j,k;
  double ****t;
  
t=(double ****) malloc((size_t)(index1*sizeof(double***)));
EXIT_IF(t == NULL,err_malloc);  
t[0]=(double ***) malloc((size_t)((index1*index2)*sizeof(double**)));
EXIT_IF(t[0] == NULL,err_malloc);  
t[0][0]=(double **) malloc((size_t)((index1*index2*index3)*sizeof(double*)));
EXIT_IF(t[0][0] == NULL,err_malloc);  
t[0][0][0]=(double *) malloc((size_t)((index1*index2*index3*index4)*sizeof(double)));
EXIT_IF(t[0][0][0] == NULL,err_malloc);  
  
 for(k=1;k<index3;k++)
   t[0][0][k]=t[0][0][k-1]+index4;
 for(i=1;i<index1;i++)
   {
     t[i]=t[i-1]+index2;
     t[i][0]=t[i-1][0]+index2*index3;
     for(j=1;j<index2;j++)
       {
	 t[i][j]=t[i][j-1]+index3;
	 t[i][j][0]=t[i][j-1][0]+index2*index3*index4;
	 for(k=1;k<index3;k++)
	   t[i][j][k]=t[i][j][k-1]+index4;
       }
   }
 
 return t;
}

void
free_dtensor4(double ****m)
{
  free((FREE_ARG) (m[0][0][0]));
  free((FREE_ARG) (m[0][0]));
  free((FREE_ARG) (m[0]));
  free((FREE_ARG) (m));
}

int ****
itensor4(long index1, long index2, long index3, long index4)
     /* allocate a int 4tensor */
{
  long i,j,k;
  int ****t;
  
t=(int ****) malloc((size_t)(index1*sizeof(int***)));
EXIT_IF(t == NULL,err_malloc);  
t[0]=(int ***) malloc((size_t)((index1*index2)*sizeof(int**)));
EXIT_IF(t[0] == NULL,err_malloc);  
t[0][0]=(int **) malloc((size_t)((index1*index2*index3)*sizeof(int*)));
EXIT_IF(t[0][0] == NULL,err_malloc);  
t[0][0][0]=(int *) malloc((size_t)((index1*index2*index3*index4)*sizeof(int)));
EXIT_IF(t[0][0][0] == NULL,err_malloc);  
  
 for(k=1;k<index3;k++)
   t[0][0][k]=t[0][0][k-1]+index4;
 for(i=1;i<index1;i++)
   {
     t[i]=t[i-1]+index2;
     t[i][0]=t[i-1][0]+index2*index3;
     for(j=1;j<index2;j++)
       {
	 t[i][j]=t[i][j-1]+index3;
	 t[i][j][0]=t[i][j-1][0]+index2*index3*index4;
	 for(k=1;k<index3;k++)
	   t[i][j][k]=t[i][j][k-1]+index4;
       }
   }
 
 return t;
}

void
free_itensor4(int ****m)
{
  free((FREE_ARG) (m[0][0][0]));
  free((FREE_ARG) (m[0][0]));
  free((FREE_ARG) (m[0]));
  free((FREE_ARG) (m));
}

#ifdef FFTW_IN_USE
fftw_complex **
fftwmatrix(long nrow, long ncol)
/* allocate a double matrix with subscript range m[0..nrow-1][0..ncol-1] */
{
  long i;
  fftw_complex **m;
  
  /* allocate pointers to rows */
  m=(fftw_complex **) fftw_malloc((size_t) (nrow*sizeof(fftw_complex*)));
  EXIT_IF(m == NULL,err_malloc);  

  /* allocate rows and set pointers to them */
m[0]=(fftw_complex *) fftw_malloc((size_t) ((nrow*ncol)*sizeof(fftw_complex)));
  EXIT_IF(m[0] == NULL,err_malloc);  
  
  for(i=1;i<nrow;i++)
    m[i]=m[i-1]+ncol;
  
  /* return pointer to array of pointers to rows */
  return m;
}

void
free_fftwmatrix(fftw_complex **m)
     /* free a double matrix allocated by dmatrix() */
{
  fftw_free((FREE_ARG) (m[0]));
  fftw_free((FREE_ARG) (m));
}

fftw_complex ***
fftwtensor3(long index1, long index2, long index3)
     /* allocate a int 3tensor */
{
  long i,j;
  fftw_complex ***t;
  
  /* allocate pointers to pointers to rows */
  t=(fftw_complex ***) fftw_malloc(index1*sizeof(fftw_complex**));
  EXIT_IF(t == NULL,err_malloc);  
  
  /* allocate pointers to rows and set pointers to them */
  t[0]=(fftw_complex **) fftw_malloc((index1*index2)*sizeof(fftw_complex*));
  EXIT_IF(t[0] == NULL,err_malloc);  
  
  /* allocate rows and set pointers to them */
t[0][0]=(fftw_complex *) fftw_malloc((index1*index2*index3)*sizeof(fftw_complex));
  EXIT_IF(t[0][0] == NULL,err_malloc);  
  
  for(j=1;j<index2;j++)
    t[0][j]=t[0][j-1]+index3;
  for(i=1;i<index1;i++)
    {
      t[i]=t[i-1]+index2;
      t[i][0]=t[i-1][0]+index2*index3;
      for(j=1;j<index2;j++)
	t[i][j]=t[i][j-1]+index3;
    }

  /* return pointer to array of pointers to rows */
  return t;
}

void
free_fftwtensor3(fftw_complex ***m)
{
  fftw_free((FREE_ARG) (m[0][0]));
  fftw_free((FREE_ARG) (m[0]));
  fftw_free((FREE_ARG) (m));
}

#endif /* FFTW_IN_USE */
