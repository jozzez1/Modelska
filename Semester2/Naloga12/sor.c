#include "sor.h"

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
//    print_array (zeta, N);
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
    print_array (psi, N);
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
//        printf ("%lf\n", *norm2);
    } while (*norm2 > p2);
}
