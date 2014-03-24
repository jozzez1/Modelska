#include "sor.h"
#include "navier_stokes.h"

void
initialize_variables (double * psi, double * zeta, double * u, double * v,
       const unsigned int N)
{
    unsigned int i,j,k;
    for (i = N; i--;)
    {
        for (j = N; j--;)
        {
            k = j + i*N;
            psi [k]   = 0;
            zeta [k]  = 0;
            u [k]     = 0;
            v [k]     = 0;
        }

        // the lowest edge for u is a bit different
        u [i + (N-1)*N] = 1;
    }
} 

void
get_vxy (double * u, double * v,
        const double * psi, const unsigned int N)
{
    unsigned int i,j,k;
    for (i = N-1; i--;)
    {
        for (j = N-1; j--;)
        {
            k = j+1 + (i+1)*N;
            u [k] = (psi[k + N] - psi[k - N]) * 0.5 * N;
            v [k] = (psi[k - 1] - psi[k + 1]) * 0.5 * N;
        }
    }
}

void
iterate_zeta (double * zeta,
        const double * u, const double * v, const double * delta,  const double * Re, const unsigned int N)
{
    unsigned int i,j,k;
    for (i = N-1; i--;)
    {
        for (j = N-1; j--;)
        {
            k = j+1 + (i+1)*N;
            zeta[k] += *delta * N * (N * (zeta[k+1] + zeta[k-1] + zeta[k+N] + zeta[k-N] - 4*zeta[k])/ *Re
                    - 0.5*(u[k+1]*zeta[k+1] - u[k-1]*zeta[k-1] + v[k+N]*zeta[k+N] - v[k-N]*zeta[k-N]));
        }
    }
}

void
fix_zeta_boundaries (double * zeta, const double * psi, const unsigned int N)
{
    unsigned int k;
    for (k = N; k--;)
    {
        // the upper boundary: zeta_ij, i = 0, j = 0,N-1
        zeta [k] = 2*pow2(N) * psi[k + N];

        // left boundary
        zeta [k*N] = 2*pow2(N) * psi[1 + k*N];

        // right boundary
        zeta [N-1 + k*N] = 2*pow2(N) * psi[N-2 + k*N];

        // bottom boundary, yay! velocity is non-zero
        zeta [k + (N-1)*N] = 2*N * (N * psi[k + (N-2)*N] - 1);
    }
}
