#include <stdio.h>
#include <unistd.h>
#include "zero_solve.h"
#include "star_system.h"

int main (int argc, char ** argv)
{
    double M1   = 2000,
           M2   = 1000,
           eps  = 0.4,
           dt   = 1e-3,
           top  = 3,
           t    = 0,
           T    = 0;

    int arg;

    while ((arg = getopt(argc, argv, "M:m:e:T:d:h")) != -1)
    {
        switch (arg)
        {
            case 'M':   M1   = atof(optarg); break;
            case 'm':   M2   = atof(optarg); break;
            case 'e':   eps  = atof(optarg); break;
            case 'T':   top  = atof(optarg); break;
            case 'd':   dt   = atof(optarg); break;
            case 'h':
                        fprintf (stdout, "-M:   mass of the 1st star\n");
                        fprintf (stdout, "-m:   mass of the 2nd star\n");
                        fprintf (stdout, "-e:   eccentricity of the star orbits\n");
                        fprintf (stdout, "-T:   how many star cycles we will compute\n");
                        fprintf (stdout, "-d:   time step duration\n");
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
