#include <stdio.h>
#include <unistd.h>
#include "zero_solve.h"
#include "star_system.h"

int main (int argc, char ** argv)
{
    double M1   = 2000,
           M2   = 1000,
           pfi  = sqrt(M1*M2)*(1 + sqrt(1 - sqrt((M1*M1 + M2*M2)/(M1*M1 + 2*M2*M1 + M2*M2)))),
           eps  = 0.4,
           dt   = 1e-3,
           top  = 3,
           t    = 0,
           T    = 0;

    int arg;

    while ((arg = getopt(argc, argv, "M:m:l:e:h")) != -1)
    {
        switch (arg)
        {
            case 'M':   M1   = atof(optarg); break;
            case 'm':   M2   = atof(optarg); break;
            case 'p':   pfi  = atof(optarg); break;
            case 'l':   eps  = atof(optarg); break;
            case 'h':
                        fprintf (stdout, "-M:   mass of the 1st star\n");
                        fprintf (stdout, "-m:   mass of the 2nd star\n");
                        fprintf (stdout, "-l:   total angular momentum of the star system\n");
                        fprintf (stdout, "-e:   eccentricity of the star orbits\n");
                        fprintf (stdout, "-h:   print this list\n");
                        exit (EXIT_SUCCESS);
            default:
                        fprintf (stderr, "Wrong program usage! Try\n%s -h\n", argv[0]);
                        exit (EXIT_FAILURE);
        }
    }

    binary sys;
    init_system (&sys, M1, M2, eps);
    T = solar_year (sys);
    top *= T;
    dt  *= T;

    double rho, phi;
    while (t < top)
    {
        trajectory (&rho, &phi, sys, t);
        printf ("%lf\t%lf\t%lf\n", t, rho, phi);
        t += dt;
    }
    return 0;
}
