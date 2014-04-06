#ifndef POINCARE_H
#define POINCARE_H
#include "integrator.h"
#include <stdio.h>

static const size_t Limit = 100;

// poincare section
void poincare (void (*scheme) (planet *, binary, params, double),
        planet * omikron, binary * sys, double dt, double precision, FILE * fout);


#endif
