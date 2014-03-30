#ifndef GLOBAL_H
#define GLOBAL_H

#include <cmath>

// a point in a 3d space
typedef struct point3
{
    double x,
           y, 
           z;
} point3;

// summation with modulo 3
unsigned int dotplus (unsigned int m, unsigned int n);

// power of 2 (e.g. square)
double pow2 (double a);

// distance between two points (Euclidian norm of the difference vector)
double dist2 (point3 a, point3 b);

#endif
