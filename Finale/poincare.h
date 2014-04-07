#ifndef POINCARE_H
#define POINCARE_H
#include "integrator.h"
#include <stdio.h>

// poincare section
void poincare (void (*scheme) (planet *, binary, params, double),
        planet * omikron, binary * sys, double dt, double precision, size_t Limit, FILE * fout);


#endif
