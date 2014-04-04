#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#define X0   1.35120719195966
#define X1  -1.70241438391932

#include <math.h>
#include "star_system.h"

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

#endif
