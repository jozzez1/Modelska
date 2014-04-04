#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#include <math.h>
#include <stdio.h>
#include "star_system.h"

// we throw in a few constants
static const double X0 =  1.3512071919596577718181151794851757586002349853515625; // X0 = 1.0 / (2 - cbrt(2))
static const double X1 = -1.7024143839193153215916254339390434324741363525390625; // X1 = cbrt(2) / (cbrt(2) - 1)

// solutions A of w for integrator S8
static const double w [8] =
{
    0,
    -0.161582374150097e+1,
    -0.244699182370524e+1,
    -0.716989419708120e-2,
    +0.244002732616735e+1,
    +0.157739928123617e+0,
    +0.182020630970714e+1,
    +0.104242620869991e+1
};

// we save our parameters in this here struct
typedef struct params
{
    double c [16],
           d [15];
} params;


// planet parameters go here
typedef struct planet
{
    double zeta,
           psi,
           p_psi,
           p_zeta;
} planet;

// I'll call the planet omikron, pun intended :P
// A_step is: x + c*dt*{x, T_zeta}/2 = ...
void A_step (planet * omikron, double dt, double c);

// B_step is: x + c*dt*{x, V_eff} = ...
void B_step (planet * omikron, binary sys, double dt, double c);

// S2 integrator consists of multiple steps
void S2 (planet * omikron, binary sys, double dt, double c);

// and finally, the S4 scheme
void S4 (planet * omikron, binary sys, double dt);

void S4a (planet * omikron, binary sys, params p, double dt);

// I'll use this one, since S4 isn't precise enough for long times
void S8 (planet * omikron, binary sys, params p, double dt);

// adaptive step control for either method, based on Richardson's extrapolation
void adaptive_step (void (* scheme) (planet *, binary, params, double),
        planet * omikron, binary * sys, params p, double dt, double * t, double precision);

// we read wA and put initialize parameters for S8
void init_params (params * p);

// solver for the orbits
void solver (planet * omikron, binary * sys, double dt, double T, FILE * fout);

// improved S8 solver
void solver_S8 (planet * omikron, binary * sys, double dt, double T, FILE * fout);

// adaptive solver
void adaptive_solver (void (* scheme) (planet *, binary, params, double),
        planet * omikron, binary * sys, double dt, double T, double precision, FILE * fout);

#endif

