#ifndef POINCARE_H
#define POINCARE_H
#include "integrator.h"
#include <stdio.h>

// poincare section
void poincare (scheme_fp scheme, planet * omikron, binary * sys, double * t_start, double dt, double precision, size_t Limit, FILE * fout);

// count the number of lines in a file
size_t linecount (FILE * fin);

// keep going from where we left off
void Continue (scheme_fp scheme, planet * omikron, binary * sys, double dt, double precision, size_t Limit, FILE * fio);

#endif
