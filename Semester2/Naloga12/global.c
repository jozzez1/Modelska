#include "global.h"

inline double
pow2 (const double x)
{
    return x*x;
}

inline void
print_array (const double * xi, unsigned long N)
{
    unsigned int i,j;
    for (i = 0; i <= N-1; i++)
    {
        for (j = 0; j <= N-1; j++)
            printf ("% .1lf\t", xi[j + i*N]);
        printf ("\n");
    }
    printf ("\n");
}

void
swap (double ** A, double ** B, const unsigned int N)
{
    double * swap = malloc (N * sizeof(double));
    swap = *B;
    *B = *A;
    A = &swap;
    free (swap);
}

