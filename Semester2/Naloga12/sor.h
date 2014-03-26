#ifndef SOR_H
#define SOR_H
#include "global.h"
#include <assert.h>
#include <math.h>

// the even SOR step
void step_even (double * psi, double * norm2, double * xi, const double * w, const double * zeta, const unsigned int N);

// the odd SOR step
void step_odd (double * psi, double * norm2, double * xi, const double * w, const double * zeta, const unsigned int N);

// next overrelaxation coefficient
void next_w (double * w, const double * J);

// first full SOR step is different, because of w_0 = 1 and w_1/2 = ...
void first_step (double * psi, double * w, double * xi, double * norm2, const double * zeta, const double * J, const unsigned int N);

// the rest of the steps are regular, so ... YAY!
void full_step (double * psi, double * w, double * xi, double * norm2, const double * zeta, const double * J, const unsigned int N);

// solve the Poisson equation
void SOR (double * psi, double * w, double * xi, double * norm2, const double * zeta, const double * J, const double p, const unsigned int N);

#endif
