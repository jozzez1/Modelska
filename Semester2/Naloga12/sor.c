#include "sor.h"

void
get_xi (double * xi, double * norm2,
        const double * psi, const double * zeta, const unsigned int N)
{
    *norm2 = 0;
    unsigned int i,j, k;
    for (i = N-2; i--;) // xi must not change the boundary!
    {
        for (j = N-2; j--;)
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

void
next_w (double * w, const double * J)
{
    *w = 1.0/(1 - 0.25*pow2(*J) * *w);
}

void
first_step (double * psi, double * w,
        const double * J, const double * xi, const unsigned int N)
{
    *w = 1;
    step_even (psi, w, xi, N);
    *w = 1.0 / (1 - 0.5*pow2(*J));
    step_odd (psi, w, xi, N);
}

void
full_step (double * psi, double * w,
        const double * J, const double * xi, const unsigned int N)
{
    next_w (w, J);
    step_even (psi, w, xi, N);
    next_w (w, J);
    step_odd (psi, w, xi, N);
}

void
SOR (double * psi, const double * zeta,
        double * xi, double * w, const double * J, double * norm2, const double p, const unsigned int N)
{
    get_xi (xi, norm2, psi, zeta, N);
    first_step (psi, w, J, xi, N);
    double p2 = pow2(p);
    do
    {
        get_xi (xi, norm2, psi, zeta, N);
        full_step (psi, w, J, xi, N);
    } while (*norm2 > p2);
}

