#include "global.h"

inline double
pow2 (const double x)
{
    return x*x;
}

inline double
mabs (const double *x)
{
    double re = *x;
    if (x < 0) re = -re;
    return re;
}

inline double
get_greater (const double *x, const double *y)
{
    double re = mabs(x),
           fe = mabs(y);
    if (re < fe) re = fe;
    return re;
}

inline void
print_array (FILE * fout, const double * xi, const unsigned long N)
{
    unsigned int i,j;
    for (i = 0; i <= N-1; i++)
    {
        for (j = 0; j <= N-1; j++)
            fprintf (fout, "% .1lf\t", xi[j + i*N]);
        fprintf (fout, "\n");
    }
    fprintf (fout, "\n");
}

void
swap (double ** A, double ** B, double ** swp)
{
    // we just swap the addresses of the arrays
    *swp = *B;
    *B = *A;
    *A = *swp;
}

