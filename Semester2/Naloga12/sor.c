#include "sor.h"

double
pow2 (const double x)
{
    return x*x;
}

void
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
get_xi (double * xi, double * norm2,
        const double * psi, const double * zeta, const unsigned int N)
{
    *norm2 = 0;
    unsigned int i,j, k;
    double h = 1.0/N;
    for (i = N-2; i--;) // xi must not change the boundary!
    {
        for (j = N-2; j--;)
        {
            k = j+1 + (i+1)*N;
            xi[k] = psi[k+1] + psi[k-1] + psi[k+N] + psi[k-N] - 4*psi[k] - zeta[k]*h*h;
            *norm2 += pow2(xi[k]);
        }
    }
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
    } while (k+2);
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
    } while (k+1);
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
    first_step (psi, w, xi, norm2, zeta, J, N);
    double p2 = pow2(p);
    do
    {
        full_step (psi, w, xi, norm2, J, xi, N);
    } while (*norm2 > p2);
}

