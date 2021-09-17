#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <float.h>

#include "pcg_basic.h"
#include "tools.h"
#include "exit_if.h"
#include "erreurs.h"
#include "allocation.h"

double maxdouble32;
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
