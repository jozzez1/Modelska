#ifndef NAVIER_STOKES_H
#define NAVIER_STOKES_H
#include "sor.h"

// initialize the arrays to their initial states
void initial_conditions (double * zeta, double * u, const unsigned int N);

// get the velocities from psi
void get_vxy (double * u, double * v, const double * psi, const unsigned int N);

// iterate zeta by one time step
void iterate_zeta (double * zeta, const double * u, const double *v, const double * delta, const double * Re, const unsigned int N);

// fix zeta boundaries
void fix_zeta_boundaries (double * zeta, const double * psi, const unsigned int N);

#endif
