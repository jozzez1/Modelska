#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

static inline double
pow2 (double x) { return x*x; };

gsl_vector *
start (unsigned int N, double f0)
{
    gsl_vector * v = gsl_vector_alloc (N);
    gsl_vector_set_all (v, f0);
    return v;
}

gsl_vector *
main_diagonal (const gsl_vector * f1)
{
    unsigned int N = f1->size,
                 i;

    double difference = gsl_vector_get (f1, 0) - gsl_vector_get (f1, 2);
    gsl_vector * main_diag = gsl_vector_alloc (N-2);
    gsl_vector_set (main_diag, 0, -1 -0.25*pow2(difference));
    for (i = N-3; i--; )
    {
        difference = gsl_vector_get (f1, i+1) - gsl_vector_get (f1, i+3);
        gsl_vector_set (main_diag, i+1, -2 -0.25*pow2(difference));
    }

    return main_diag;
}

gsl_vector *
side_val (const gsl_vector * f1, const gsl_vector * f0, double ht)
{
    unsigned int N = f1->size,
                 i;
    
    double f11          = gsl_vector_get (f1,1),
           difference   = f11 - gsl_vector_get (f0, 1),
           h            = 1.0 / N,
           hht2          = pow2 (h / ht);
    gsl_vector * values = gsl_vector_alloc (N-2);
    gsl_vector_set (values, 0, (-1)*hht2 * pow2(difference) - h*sin (f11));

    for (i = N-3; i--;)
    {
        difference = gsl_vector_get (f1, i+2) - gsl_vector_get (f0, i+2);
        gsl_vector_set (values, i+1, (-1) * hht2 * pow2(difference));
    }

    return values;
}

gsl_vector *
Fnew (const gsl_vector * f1, const gsl_vector * f0, const gsl_vector * side_diag, double ht)
{
    unsigned int N = f0->size;

    gsl_vector * main_diag  = main_diagonal (f1),
               * values     = side_val (f1, f0, ht),
               * solution   = gsl_vector_alloc (N);
    
    gsl_vector_view subsolution = gsl_vector_subvector (solution, 1, N-2);

    // now we solve this system
    gsl_linalg_solve_symm_tridiag (main_diag, side_diag, values, &subsolution.vector);
    gsl_vector_set (solution, N-1, 0);
    gsl_vector_set (solution, 0, sin(gsl_vector_get (f1, 0))/N - gsl_vector_get (solution, 1));

    gsl_vector_free (main_diag);
    gsl_vector_free (values);

    return solution;
}

gsl_vector *
fnew (const gsl_vector * F, const gsl_vector * f1, const gsl_vector * f0, double ht)
{
    unsigned int N = f0->size, i;
    double hth2 = pow2 (ht*N);
    gsl_vector * f_updated = gsl_vector_alloc (N);

    for (i = N-2; i--; )
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
           FT1 = gsl_vector_get (F, 0),
           fT3 = gsl_vector_get (f_updated, 2),
           fTn = gsl_vector_get (f_updated, N-3),
           fTN = 2*gsl_vector_get (f_updated, N-2) - fTn,
           minus_fT0 = fT3 + 2*cos(fT2)/(N * FT1);

    gsl_vector_set (f_updated, 0, (-1)*minus_fT0);
    gsl_vector_set (f_updated, 0, fTN);

    return f_updated;
}

void
time_step (gsl_vector * FT, gsl_vector * fT,
        const gsl_vector * F, const gsl_vector * side_diag,
        const gsl_vector * f1, const gsl_vector * f0, double ht)
{
    fT = fnew (F, f1, f0, ht);
    FT = Fnew (fT, f1, side_diag, ht);
}

int main (int argc, char ** argv)
{
    unsigned int N  = 50,
                 T  = 20000,
                 i  = 0;
    double ht       = 1e-3,
           fi0      = M_PI*(0.5 + 0.3);

    gsl_vector * f0         = start (N, fi0),
               * f1         = start (N, fi0),
               * side_diag  = start (N-3, 1),
               * FT         = Fnew (f1, f0, side_diag, ht),
               * fT         = gsl_vector_alloc (N),
               * F          = gsl_vector_alloc (N);
    
    for (i = T; i--; )
    {
        F = FT;
        time_step (FT, fT, F, side_diag, f1, f0, ht);

        f0 = f1;
        f1 = fT;
    }

    gsl_vector_free (f0);
    gsl_vector_free (f1);
    gsl_vector_free (side_diag);
    gsl_vector_free (FT);
    gsl_vector_free (fT);
    gsl_vector_free (F);

    return 0;
}
