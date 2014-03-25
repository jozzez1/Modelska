#ifndef GLOBAL_H
#define GLOBAL_H
#include <stdio.h>
#include <stdlib.h>

// power of 2
double pow2 (const double x);

// for debugging
void print_array (const double * xi, unsigned long N);

// swap the two arrays in constant time
void swap (double ** A, double ** B, const unsigned int N);

#endif
