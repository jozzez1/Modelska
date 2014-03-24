#include "sor.h"

void
get_xi (double * xi, double * norm2,
        const double * psi, const double * zeta, const unsigned int N)
{
    *norm2 = 0;
    unsigned int i,j, k;
    for (i = N-1; i--;)
    {
        for (j = N-1; j--;)
        {
            k = j+1 + (i+1)*N;
            xi[k] = psi[k+1] + psi[k-1] + psi[k+N] + psi[k-N] - 4*psi[k] - zeta[k];
            *norm2 += pow2(xi[k]);
        }
    }
}

void
step_even (double * psi,
        const double * w, const double * xi, const unsigned int N)
{
    unsigned int i,j,k;
    for (i = N; i -= 2;)
    {
        for (j = N; j -= 2;)
        {
            k = j + i*N;
            psi[k] += 0.25 * (*w) * xi[k];
        }
    }
}

void
step_odd (double * psi,
        const double * w, const double * xi, const unsigned int N)
{
    unsigned int i,j,k;
    for (i = N; i -= 2;)
    {
        for (j = N; j -= 2;)
        {
            k = j+1 + (i+1)*N;
            psi[k] += 0.25 * (*w) * xi[k];
        }
    }
}
