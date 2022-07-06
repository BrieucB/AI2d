#include "pcg_basic.h"
#ifndef TOOLS_H
#define TOOLS_H

#define M_PI 3.14159265358979323846

int WhichBox(int, double, double);

void init_frame(int, int, int *);

void shuffleinbox2d(int, int, int,
	       double *, double *,
		    int *, int *, int *, int *);


void computeDisplacement(int , int , int , double , double , double , double , double ,
			 double , pcg32_random_t *, double *, double *, int *, int *,
			 int *, int *, int *, double *, double *);

int updatePositions(int , int , int , double *, double *, double *, double *,
		    int *, int *);

void constant(void);

double genrand32_real1(pcg32_random_t *);

double genrand32_real2(pcg32_random_t *);

void BoxMuller(pcg32_random_t *, double *, double *);

#endif
