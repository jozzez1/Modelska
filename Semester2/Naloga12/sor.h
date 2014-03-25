#ifndef SOR_H
#define SOR_H

static inline double
pow2 (double x) { return x*x; };

// get XI and its norm
void get_xi (double * xi, double * norm2, const double * psi, const double * zeta, const unsigned int N);

// the even SOR step
void step_even (double * psi, const double * w, const double * xi, unsigned int N);

// the odd SOR step
void step_odd (double * psi, const double * w, const double * xi, unsigned int N);

// next overrelaxation coefficient
void next_w (double * w, const double * J);

// first full SOR step is different, because of w_0 = 1 and w_1/2 = ...
void first_step (double * psi, double * w, double * xi, double * norm2, const double * zeta, const double * J, const unsigned int N);

// the rest of the steps are regular, so ... YAY!
void full_step (double * psi, double * w, double * xi, double * norm2, const double * zeta, const double * J, const unsigned int N);

// solve the Poisson equation
void SOR (double * psi, double * w, double * xi, double * norm2, const double * zeta, const double * J, const double p, const unsigned int N);

// for debugging purposes
#include <stdio.h>
void print_array (const double * xi, unsigned long N);

#endif
