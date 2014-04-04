#include <stdio.h>
#include <unistd.h>
#include "zero_solve.h"
#include "star_system.h"
#include "integrator.h"

int main (int argc, char ** argv)
{
    double M1   = 2000,
           M2   = 1000,
           eps  = 0.4,
           dt   = 1e-3,
           top  = 3,
           T    = 0,
           zeta = 20,
           psi  = 0,
           p_zeta = 0,
           p_psi  = 100,
           prec   = 1e-3;

    int option = 0,
        arg;

    while ((arg = getopt(argc, argv, "M:m:e:T:d:z:l:p:r:O:h")) != -1)
    {
        switch (arg)
        {
            case 'M':   M1   = atof(optarg); break;
            case 'm':   M2   = atof(optarg); break;
            case 'e':   eps  = atof(optarg); break;
            case 'T':   top  = atof(optarg); break;
            case 'd':   dt   = atof(optarg); break;
            case 'z':   zeta = atof(optarg); break;
            case 'r':   prec = atof(optarg); break;
            case 'l':   p_psi   = atof(optarg); break;
            case 'p':   p_zeta  = atof(optarg); break;
            case 'O':   option  = atoi(optarg); break;
            case 'h':
                        fprintf (stdout, "-M:   mass of the 1st star\n");
                        fprintf (stdout, "-m:   mass of the 2nd star\n");
                        fprintf (stdout, "-e:   eccentricity of the star orbits\n");
                        fprintf (stdout, "-T:   how many star cycles we will compute\n");
                        fprintf (stdout, "-d:   time step duration\n");
                        fprintf (stdout, "-z:   zeta of the planet\n");
                        fprintf (stdout, "-l:   starting angular momentum of the planet\n");
                        fprintf (stdout, "-p:   starting radial momentum of the planet\n");
                        fprintf (stdout, "-r:   precision for the adaptive step\n");
                        fprintf (stdout, "-o:   option -- which integrator to choose\n");
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

    planet omikron =
    {
        .zeta = zeta,
        .psi    = psi,
        .p_psi  = p_psi,
        .p_zeta = p_zeta
    };

    switch (option)
    {
        case 0: solver (&omikron, &sys, dt, top, stdout);    break;
        case 1: solver_S8 (&omikron, &sys, dt, top, stdout); break;
        case 2: adaptive_solver (&S8, &omikron, &sys, dt, top, prec, stdout); break;
        default:
                fprintf (stderr, "Unknown option\n");
                exit (EXIT_FAILURE);
    }

    return 0;
}
