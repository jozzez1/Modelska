#ifndef NAVIER_STOKES_H
#define NAVIER_STOKES_H
#include "global.h"

// initialize the arrays to their initial states
void initial_conditions (double * zeta, double * u, const unsigned int N);

// get the velocities from psi
void get_vxy (double * u, double * v, const double * psi, const unsigned int N);

// get the velocities and anso predict delta
void get_vxy_and_delta (double * u, double * v, double * delta, const double * psi, const unsigned int N);

// iterate zeta by one time step
void iterate_zeta (double * zeta, const double * tmp_zeta, const double * u, const double *v, const double * delta, const double * Re, const unsigned int N);

// fix zeta boundaries
void fix_zeta_boundaries (double * zeta, const double * psi, const unsigned int N);

#endif
