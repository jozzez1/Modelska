#ifndef STAR_SYSTEM_H
#define STAR_SYSTEM_H
// we look for different
#include <assert.h>
#include <math.h>
#include "zero_solve.h"

// struct with all the parameters of the system
typedef struct binary
{
    double epsilon, // eccentricity of the orbits
           alpha,   // potential well strength
           mu,      // reduced mass of the system
           M1,      // mass of the 1st star
           M2,      // mass of the 2nd star
           p_phi;   // angular momentum
} binary;

// initialize the system
void init_system (binary * sys, double M1, double M2, double p_phi, double p_rho);

// calculate stuff from parameters and time, like the trajectory :D
void trajectory (double * rho, double * phi, binary sys, double t);

// return the period of such a binary system
double solar_year (binary sys);

#endif
