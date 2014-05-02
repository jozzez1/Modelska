#ifndef POINCARE_H
#define POINCARE_H
#include "integrator.h"
#include <stdio.h>

#define PLANET_ESCAPE 42

// poincare section
int poincare (scheme_fp scheme, planet * omikron, binary * sys, double * t_start, double dt, double precision, size_t Limit, FILE * fout);

// count the number of lines in a file
size_t linecount (FILE * fin);

// keep going from where we left off
void Continue (scheme_fp scheme, planet * omikron, binary * sys, double dt, double precision, size_t Limit, FILE * fio);

// use bisection to find the equilibrium point
double bisect (scheme_fp scheme, planet * omikron, binary * sys, double Z, double dt, double precision);

// we can use this function to scan over a certain region in angular momentum to find equilibrium points
void equilibrium (scheme_fp scheme, planet * omikron, binary * sys, double Lmax, double dt, double precision, FILE * fout);

#endif
