#include "sor.h"

void
step_even (double * psi, double * norm2, double * xi,
        const double * w, const double * zeta, const unsigned int N)
{
    double N2 = N*N;
    unsigned int i,j,k;
    i = N-2;
    while (i+1)
    {
        j = N-2;
        while (j+1)
        {
            k = j + i*N;
            *xi = psi[k+1] + psi[k-1] + psi[k+N] + psi[k-N] - 4*psi[k] - zeta[k]/N2;
            *norm2 += (*xi)*(*xi);
            psi[k] += 0.25 * (*w) * (*xi);
            j -= 2;
        }
        i -= 2;
    }

    i = N-3;
    while (i)
    {
        j = N-3;
        while (j)
        {
            k = j + i*N;
            *xi = psi[k+1] + psi[k-1] + psi[k+N] + psi[k-N] - 4*psi[k] - zeta[k]/N2;
            *norm2 += (*xi)*(*xi);
            psi[k] += 0.25 * (*w) * (*xi);
            j -= 2;
        }
        i -= 2;
    }
    assert (isnormal(*norm2));
}

void
step_odd (double * psi, double * norm2, double * xi,
        const double * w, const double * zeta, const unsigned int N)
{
    double N2 = N*N;
    unsigned int i,j,k;
    i = N-2;
    while (i+1)
    {
        j = N-3;
        while (j)
        {
            k = j + i*N;
            *xi = psi[k+1] + psi[k-1] + psi[k+N] + psi[k-N] - 4*psi[k] - zeta[k]/N2;
            *norm2 += (*xi)*(*xi);
            psi[k] += 0.25 * (*w) * (*xi);
            j -= 2;
        }
        i -= 2;
    }

    i = N-3;
    while (i)
    {
        j = N-2;
        while (j+1)
        {
            k = j + i*N;
            *xi = psi[k+1] + psi[k-1] + psi[k+N] + psi[k-N] - 4*psi[k] - zeta[k]/N2;
            *norm2 += (*xi)*(*xi);
            psi[k] += 0.25 * (*w)*(*xi);
            j -= 2;
        }
        i -= 2;
    }
    assert (isnormal(*norm2));
}

void
next_w (double * w, const double * J)
{
    *w = 1.0/(1 - 0.25*(*J)*(*J)*(*w));
}

void
first_step (double * psi, double * w, double * xi, double * norm2, 
        const double * zeta, const double * J,const unsigned int N)
{
    *norm2 = 0;
    *w = 1;
    step_even (psi, norm2, xi, w, zeta, N);
    *w = 1.0 / (1 - 0.5*(*J)*(*J));
    step_odd (psi, norm2, xi, w, zeta, N);
}

void
full_step (double * psi, double * w, double * xi, double * norm2,
        const double * zeta, const double * J,const unsigned int N)
{
    *norm2 = 0;
    next_w (w, J);
    step_even (psi, norm2, xi, w, zeta, N);
    next_w (w, J);
    step_odd (psi, norm2, xi, w, zeta, N);
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

