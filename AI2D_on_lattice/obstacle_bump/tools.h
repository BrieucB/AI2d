#include "pcg_basic.h"
#ifndef TOOLS_H
#define TOOLS_H

#define M_PI 3.14159265358979323846

// double maxdouble32;
void constant(void);

double genrand32_real1(pcg32_random_t *);

double genrand32_real2(pcg32_random_t *);

double computeLocalQuantities(int , int , int , int *, int *, int *, int **, int **);

void updatePositions(int, int, int *, int *, int * , int , int **, int **,  double , double , double , double , double , pcg32_random_t *);

void updatePositionsObstacle(int , int , int *, int *, int *, int , int **, int **, double , double , double , double , double , pcg32_random_t *, int , int , double);

void updatePositionsObstacleBump(int , int , int *, int *, int *, int , int **, int **, double , double , double , double , double , pcg32_random_t *, int , int , double);

#endif
