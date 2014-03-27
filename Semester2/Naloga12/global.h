#ifndef GLOBAL_H
#define GLOBAL_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <mgl2/data.h>
#include <mgl2/mgl.h>

// power of 2 (e.g. square)
double pow2 (const double x);

// absolute value of a floating number
double mabs (const double *x);

// get the greater of two numbers
double get_greater (const double *x, const double *y);

// for debugging
void print_array (FILE * fout, const double * xi, const unsigned long N);

// swap the two arrays in constant time
void swap (double ** A, double ** B, double ** swp);

#endif
