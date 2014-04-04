#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#include <math.h>
#include <stdio.h>
#include "star_system.h"

// we throw in a few constants
static const double X0 =  1.35120719195966; // X0 = 1.0 / (2 - cbrt(2))
static const double X1 = -1.70241438391932; // X1 = cbrt(2) / (cbrt(2) - 1)

// planet parameters go here
typedef struct planet
{
    double zeta,
           psi,
           p_psi,
           p_zeta;
} planet;

// I'll call the planet omikron, pun intended :P
// A_step is: x + c*dt*{x, T}/2 = ...
void A_step (planet * omikron, double dt, double c);

// B_step is: x + c*dt*{x, V} = ...
void B_step (planet * omikron, binary sys, double dt, double c);

// S2 integrator consists of multiple steps
void S2 (planet * omikron, binary sys, double dt, double c);

// and finally, the S4 scheme
void S4 (planet * omikron, binary sys, double dt);

// solver for the orbits
void solver (planet * omikron, binary * sys, double dt, double T, FILE * fout);

#endif
