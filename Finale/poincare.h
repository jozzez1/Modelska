#ifndef POINCARE_H
#define POINCARE_H
#include "integrator.h"
#include <stdio.h>

// poincare section
void poincare (void (*scheme) (planet *, binary, params, double),
        planet * omikron, binary * sys, double dt, double precision, size_t Limit, FILE * fout);

size_t linecount (FILE * fin);

// keep going from where we left off
void Continue (void (*scheme) (planet *, binary, params, double),
        planet * omikron, binary * sys, double dt, double precision, size_t Limit, FILE * fio);

#endif
