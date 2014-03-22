#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

static inline double
pow2 (double x) { return x*x; };

void
tridiag_params (gsl_vector * main_diag, gsl_vector * values,
        gsl_vector * f1, gsl_vector * f0, double ht)
{
    unsigned int N = f1->size,
                 i;

    double f10          = gsl_vector_get (f1,0),
           differenceV  = gsl_vector_get(f1,1) - gsl_vector_get (f0, 1),
           differenceM  = gsl_vector_get (f1, 0) - gsl_vector_get (f1, 2),
           h            = 1.0 / N,
           hht2         = pow2 (h / ht);
    gsl_vector_set (values, 0, (-1)*hht2 * pow2(differenceV) - h*sin (f10));
    gsl_vector_set (main_diag, 0, -1 -0.25*pow2(differenceM));

    for (i = N-3; i--;)
    {
        differenceV = gsl_vector_get (f1, i+2) - gsl_vector_get (f0, i+2);
        differenceM = gsl_vector_get (f1, i+1) - gsl_vector_get (f1, i+3);
        gsl_vector_set (values, i+1, (-1)*hht2 * pow2(differenceV));
        gsl_vector_set (main_diag, i+1, -2 -0.25*pow2(differenceM));
    }
}

void
Fnew (gsl_vector * solution, gsl_vector * main_diag, gsl_vector * values,
        gsl_vector * f1, gsl_vector * f0, gsl_vector * side_diag, double ht)
{
    unsigned int N = f0->size;
    tridiag_params (main_diag, values, f1, f0, ht);
    gsl_vector_view subsolution = gsl_vector_subvector (solution, 1, N-2);

    // now we solve this system
    gsl_linalg_solve_symm_tridiag (main_diag, side_diag, values, &subsolution.vector);
    gsl_vector_set (solution, N-1, 0);
    gsl_vector_set (solution, 0, sin(gsl_vector_get (f1, 0))/N - gsl_vector_get (solution, 1));
}

void
fnew (gsl_vector * f_updated, gsl_vector * F, gsl_vector * f1, gsl_vector * f0, double ht)
{
    unsigned int N = f0->size, i;
    double hth2 = pow2 (ht*N);

    for (i = N-2; i--;)
    {
        double f11 = gsl_vector_get (f1, i+2),
               f1n = gsl_vector_get (f1, i+1),
               f10 = gsl_vector_get (f1, i),
               f00 = gsl_vector_get (f0, i+1),
               F1  = gsl_vector_get (F, i+2),
               F0  = gsl_vector_get (F, i),
               Fn  = gsl_vector_get (F, i+1),
               com = hth2 * (0.5*(f11 - f10)*(F1 - F0) + Fn*(f11 + f10 - 2*f1n)) + 2*f1n - f00;

        gsl_vector_set (f_updated, i+1, com);
    }
    double fT2 = gsl_vector_get (f_updated, 1),
           FT1 = gsl_vector_get (F, 1),
           fT3 = gsl_vector_get (f_updated, 2),
           fTn = gsl_vector_get (f_updated, N-3),
           fTN = 2*gsl_vector_get (f_updated, N-2) - fTn,
           minus_fT0 = 2*cos(fT2)/(N * FT1) + fT3;

    gsl_vector_set (f_updated, 0, minus_fT0);
    gsl_vector_set (f_updated, N-1, fTN);
}

void
time_step (gsl_vector * FT, gsl_vector * fT, gsl_vector * main_diag, gsl_vector * values,
        gsl_vector * F, gsl_vector * side_diag,
        gsl_vector * f1, gsl_vector * f0, double ht)
{
    fnew (fT, F, f1, f0, ht);
    Fnew (FT, main_diag, values, fT, f1, side_diag, ht);
}

int main (int argc, char ** argv)
{
    unsigned int N  = 20,
                 T  = 100000,
                 i;
    double ht       = 1e-4,
           fi0      = M_PI*(0.5 + 0.3);

    // allocation phase
    gsl_vector * f0         = gsl_vector_alloc (N),
               * f1         = gsl_vector_alloc (N),
               * FT         = gsl_vector_alloc (N),
               * fT         = gsl_vector_alloc (N),
               * F          = gsl_vector_alloc (N),
               * main_diag  = gsl_vector_alloc (N-2),
               * side_diag  = gsl_vector_alloc (N-3),
               * values     = gsl_vector_alloc (N-2);

    // now we fill them
    gsl_vector_set_all (f0, fi0);
    gsl_vector_set_all (f1, fi0);
    gsl_vector_set_all (side_diag, 1);
    Fnew (FT, main_diag, values, f1, f0, side_diag, ht);

    for (i = T; i--; )
    {
        gsl_vector_memcpy (F, FT);
        time_step (FT, fT, main_diag, values, F, side_diag, f1, f0, ht);
        gsl_vector_memcpy (f0, f1);
        gsl_vector_memcpy (f1, fT);
    }

    gsl_vector_add_constant(fT, -0.5*M_PI);
    gsl_vector_fprintf (stdout, fT, "%f");

    gsl_vector_free (f0);
    gsl_vector_free (f1);
    gsl_vector_free (FT);
    gsl_vector_free (fT);
    gsl_vector_free (F);
    gsl_vector_free (side_diag);
    gsl_vector_free (main_diag);
    gsl_vector_free (values);

    return 0;
}
