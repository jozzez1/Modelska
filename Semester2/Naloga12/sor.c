#include "sor.h"

void
get_xi (double * xi, double * norm2,
        const double * psi, const double * zeta, const unsigned int N)
{
    *norm2 = 0;
    unsigned int i,j, k;
    double N2 = N*N;
    for (i = N-2; i--;) // xi must not change the boundary!
    {
        for (j = N-2; j--;)
        {
            k = j+1 + (i+1)*N;
            xi[k] = psi[k+1] + psi[k-1] + psi[k+N] + psi[k-N] - 4*psi[k] - zeta[k]/N2;
            *norm2 += pow2(xi[k]);
        }
    }
    assert (!isinf(*norm2));
    assert (!isnan(*norm2));
}

void
step_even (double * psi,
        const double * w, const double * xi, const unsigned int N)
{
    unsigned int k = N*N-2;
    do
    {
        psi[k] += 0.25 * (*w) * xi[k];
        k -= 2;
    } while (k+1);
}

void
step_odd (double * psi,
        const double * w, const double * xi, const unsigned int N)
{
    unsigned int k = N*N-1;
    do
    {
        psi[k] += 0.25 * (*w) * xi[k];
        k -= 2;
    } while (k+2);
}

void
next_w (double * w, const double * J)
{
    *w = 1.0/(1 - 0.25*pow2(*J) * (*w));
}

void
first_step (double * psi, double * w, double * xi, double * norm2, 
        const double * zeta, const double * J,const unsigned int N)
{
    *w = 1;
    get_xi (xi, norm2, psi, zeta, N);
    step_even (psi, w, xi, N);
    *w = 1.0 / (1 - 0.5*pow2(*J));
    get_xi (xi, norm2, psi, zeta, N);
    step_odd (psi, w, xi, N);
}

void
full_step (double * psi, double * w, double * xi, double * norm2,
        const double * zeta, const double * J,const unsigned int N)
{
    next_w (w, J);
    get_xi (xi, norm2, psi, zeta, N);
    step_even (psi, w, xi, N);
    next_w (w, J);
    get_xi (xi, norm2, psi, zeta, N);
    step_odd (psi, w, xi, N);
}

void
SOR (double * psi, double * w, double * xi, double * norm2,
        const double * zeta, const double * J, const double p, const unsigned int N)
{
    double p2 = p*p;
    first_step (psi, w, xi, norm2, zeta, J, N);
    do
    {
        full_step (psi, w, xi, norm2, zeta, J, N);
    } while (p2 < *norm2);
}

