#ifndef ZERO_SOLVE_H
#define ZERO_SOLVE_H
// this header includes functions for root solving of our solar system
#include <math.h>
#include <stdlib.h>

// the modulo for doubles -- I took it from Seminarji/DinAna/hed.h
double mod (double a, double b);

// this function returns t*a^2/p^3, so not explicit time
double tap (double epsilon, double phi);

// this function returns 1st derivative of function tap
double dtap (double epsilon, double phi);

// and this one returns the 2nd derivative
double ddtap (double epsilon, double phi);

// and here we get phi for given t
void get_phi (double * phi, double epsilon, double t);

#endif
