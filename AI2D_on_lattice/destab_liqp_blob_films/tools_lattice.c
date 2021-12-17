#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <pcg_variants.h>
#include <time.h>
#include <float.h>

#include "tools.h"
/*
int nextMove(int np,
	     int nm,
	     double W,
	     double D,
	     double lambda,
	     double beta)
{
  
}
*/

/* RANDOM TOOLS */

void constant(void)
{
  	int i;
  
  	maxdouble64=1.;
  
  	for (i=0 ; i<64 ; i++)
    	maxdouble64 *= 2.;
    
  	return;
}

/* random double in [0;1] */
double genrand64_real1(pcg64_random_t *rng)
{
  	return (((double) pcg64_random_r(rng))/(maxdouble64-1.)) ;
}

/* random double in [0;1[ */
double genrand64_real2(pcg64_random_t *rng)
{
  	return (((double) pcg64_random_r(rng))/maxdouble64) ;
}

void BoxMuller(pcg64_random_t *rng, double *x, double *y)
{
  	double norm, spd2, spdx, spdy;
  
  	spd2=2.;
  	while ((spd2>1.)||(spd2<DBL_EPSILON))
    {
      	spdx=2.*genrand64_real1(rng)-1.;
      	spdy=2.*genrand64_real1(rng)-1.;

      	spd2=spdx*spdx+spdy*spdy;
    }
  	norm=sqrt(-2.*log(spd2)/spd2);
      
  	*x = spdx*norm;
  	*y = spdy*norm;

  	return ;
}

/* SORT methods */

void swapterm(int n, int m, double *table, int *indx)
{
  double temp1, temp2;
  
  temp1= table[n];
  table[n] = table[m];
  table[m] = temp1;

  temp2=indx[n];
  indx[n]=indx[m];
  indx[m]=temp2;
  
  return ;
}

double median(double a, double b, double c)/* mediane de (a,b,c) */
{
  double result;

  if (a<b)
    {
      if (b<c)
        result=b;
      else
        {
          if (a<c)
            result=c;
          else
            result=a;
        }
    }
  else
    {
      if (a<c)
        result=a;
      else
        {
          if (b<c)
            result=c;
          else
            result=b;
        }      
    }
  
  return result ;
}


void quickSort(int start, int end, double *table, int *indx)
{
  
  int left, right;
  double pivot;
  
  /* Si le table est de longueur nulle, il n'y a rien à faire. */
  /* end : taille du tableau a parcourir selon la convention du C */
  if(end-start <= 1)
    return ;
  
  pivot = median(table[start], table[end-1], table[(start+end-1)/2]);
  
  left=start-1;
  right=end;

  /* Sinon, on parcourt le tableau, une fois de droite à gauche, et
     une autre de gauche à droite, à la recherche d'éléments mal
     placés, que l'on permute. Si les deux parcours se croisent, on
     arrête. */

  for(;;)
    {
        do
          right--;
        while(table[right] > pivot);
        
        do left++;
        while(table[left] < pivot);

        if(left < right)
          swapterm(left, right, table, indx);
        else
          break;
    }
  
  /* Maintenant, tous les éléments inférieurs au pivot sont avant ceux
     supérieurs au pivot. On a donc deux groupes de cases à trier. On utilise
     pour cela... la méthode quickSort elle-même ! */
  right++;  
  quickSort(right, end, table, indx);
  
  quickSort(start, right, table, indx);
  
}

void sedgeSort(int start, int end, int limit, double *table, int *indx)
{
  int left, right, width;
  double pivot;
  
  width=end-start;
  
  /* Si le table est de longueur nulle, il n'y a rien à faire. */
  /* end : taille du tableau a parcourir selon la convention du C */
  if(width <= 1)
    return ;
  
  if (width < limit)
    {
      insertSort(start, end, table, indx);
      return ;
    }
  else
    {
      pivot = median(table[start], table[end-1], table[(start+end-1)/2]);
    
      left=start-1;
      right=end;

      /* Sinon, on parcourt le tableau, une fois de droite à gauche,
     et une autre de gauche à droite, à la recherche d'éléments mal
     placés, que l'on permute. Si les deux parcours se croisent, on
     arrête. */

      for(;;)
        {
          do
            right--;
          while(table[right] > pivot);
          
          do left++;
          while(table[left] < pivot);
          
          if(left < right)
            swapterm(left, right, table, indx);
          else
            break;
        }
    }/*end if limit */
  
  
  /* Maintenant, tous les éléments inférieurs au pivot sont avant ceux
     supérieurs au pivot. On a donc deux groupes de cases à trier. On utilise
     pour cela... la méthode sedgeSort elle-même ! */

  right++;
  
  sedgeSort(right, end, limit, table, indx);
  sedgeSort(start, right, limit, table, indx);
  
}


void insertSort(int start, int end, double *table, int *indx)
{
  int i;
  
  for (i=start+1 ; i<end ; i++)
    insert(start, i+1, table[i], indx[i], table, indx);

  return ;
}


void insert(int start, int end, double value, int rank, double *table, int *indx)
{
  int j;

  /* end : taille du tableau a parcourir selon la convention du C */

  end--;
  
  for (j = end; j > start && table[j-1] > value; j--)
    {
      table[j] = table[j-1];
      indx[j] = indx[j-1];
    }
  
    table[j] = value;
    indx[j] = rank;
    
    return;
}

/* DICHOTOMY methods */

int dichotomyStep(int imin,
		  int imax,
		  double *arr,
		  double x)
{
  int i = (imin+imax) / 2; 

  if (arr[i] < x)
    {
      imin=i;
    }
  
  else
    {
      if (arr[i-1] > x)
	{
	  imax=i;
	}
      else
	{
	  return i-1;
	}
    }
  return dichotomyStep(imin, imax, arr, x); 
} 

int dichotomySearch(double *vec,
		    int n,
		    double x)
{
  return dichotomyStep(0, n, vec, x);
}


void printIntVec(int *vec,
		 int n)
{
  int k;
  printf("[");
  for(k=0; k<n-1; k++)
    {
      printf("%d, ", vec[k]);
    }
  printf("%d]\n", vec[n-1]);
  return;
}


void printDoubleVec(double *vec,
		    int n)
{
  int k;
  printf("[");
  for(k=0; k<n-1; k++)
    {
      printf("%f, ", vec[k]);
    }
  printf("%f]\n", vec[n-1]);
  return;
}

double * cumSum(double *vec,
		int n)
{
  double *res;
  int k;
  
  res = (double *)malloc((n+1)*sizeof(double));
  res[0]=0;

  for (k=0 ; k<n ; k++)
    {
      res[k+1]=res[k]+vec[k];
    }
  return res;
}

int * dichoTowerSampling(int N,
			 double *P,
			 int n,
			 pcg64_random_t *rng)
{
  int *X;
  int i;
  double *C;
  double x;
  
  C = cumSum(P, n); // C is of size n+1, C[0]=0
  printDoubleVec(P, n);
  
  X = (int *)malloc(N*sizeof(int));

  for (i=0 ; i<N ; i++)
    {
      x=genrand64_real2(rng);
      X[i]=dichotomySearch(C, n+1, x);
    }
  return X;
}


void getMagnetization(int *vec_pos,
		      int *vec_spin,
		      int *vec_M,
		      int L,
		      int N)
{
  int k;
  
  for(k=0 ; k<L ; k++) // reset vec_M 
    {
      vec_M[k]=0;
    }
  
  for(k=0 ; k<N ; k++)
    {
      vec_M[vec_pos[k]]+=vec_spin[k];
    }

  return;
}


