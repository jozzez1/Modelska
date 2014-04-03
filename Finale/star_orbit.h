#ifndef STAR_ORBIT_H
#define STAR_ORBIT_H
// we look for different
#include <assert.h>
#include <math.h>
#include "zero_solve.h"

// struct with all the parameters of the system
typedef struct system
{
    double epsilon, // eccentricity of the orbits
           alpha,   // potential well strength
           E,       // energy
           mu,      // reduced mass of the system
           M1,      // mass of the 1st star
           M2,      // mass of the 2nd star
           p_phi,   // angular momentum
           p_rho;   // radial momentum
} system;

// initialize the system
void init_system (system * u, double M1, double M2, double p_phi, double p_rho);

// calculate stuff from parameters and time, like the trajectory :D
void trajectory (double * rho, double * phi, system sys, double t);

#endif
