#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "global.h"
#include "sor.h"
#include "navier_stokes.h"

void solve (const unsigned int N,
        const unsigned long T,
        const unsigned long frames,
        const double delta,
        const double Re,
        const double precision)
{
    // variable declarations
    // --------------------------------------------------------
    double * tmp    = (double *) calloc (N*N, sizeof(double)),
           * zeta   = (double *) calloc (N*N, sizeof(double)),
           * psi    = (double *) calloc (N*N, sizeof(double)),
           * xi     = (double *) calloc (N*N, sizeof(double)),
           * u      = (double *) calloc (N*N, sizeof(double)),
           * v      = (double *) calloc (N*N, sizeof(double));

    double J        = cos(M_PI/N),
           norm2    = 0,
           w        = 0;

    // solution
    // --------------------------------------------------------

    // first we set the non-zero conditions
    initial_conditions (zeta, u, N);
    // and now we iterate
//    unsigned int i;
//    for (i = frames; i--;)
//    {
        SOR (psi, &w, xi, &norm2, zeta, &J, precision, N);
        get_vxy (u, v, psi, N);
        print_array (tmp, N);
        swap (&zeta, &tmp, N*N);
        print_array (tmp, N);
        iterate_zeta (zeta, tmp, u, v, &delta, &Re, N);
        print_array (zeta, N);
//        fix_zeta_boundaries (zeta, psi, N);
//    }

    // variable deallocation
    // --------------------------------------------------------
    free (zeta);
    free (psi);
    free (xi);
    free (u);
    free (v);
}

int main (int argc, char ** argv)
{
    unsigned int N  = 11;
    unsigned long T = 1000,
                  F = 100;
    double delta    = 1e-3,
           Re       = 1e-2,
           precision= 1e-5;

    assert (N & 1);     // SOR converges only if N is an odd number
    assert (delta < 0.4/N);

    solve (N, T, F, delta, Re, precision);

    return 0;
}
