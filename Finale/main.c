#include <stdio.h>
#include "zero_solve.h"

int main (void)
{
    double eps = 0.4,
           t = 0,
           phi = 0;

    while (t < 100)
    {
        get_phi (&phi, eps, t);
        printf ("%lf\t%lf\n", t, phi);
        t += 0.003;
    }
    return 0;
}
