#include <stdio.h>
#include <unistd.h>
#include "zero_solve.h"
#include "star_system.h"
#include "integrator.h"
#include "poincare.h"

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
        out    = 0,
        arg;

    FILE * fout;
    char * filename = NULL;
    while ((arg = getopt(argc, argv, "M:m:e:T:d:z:l:p:r:O:o:h")) != -1)
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
            case 'o':   
                        filename = optarg;
                        out      = 1;           break;
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
                        fprintf (stdout, "-O:   option -- which integrator to choose\n");
                        fprintf (stdout, "-h:   print this list\n");
                        exit (EXIT_SUCCESS);
            default:
                        fprintf (stderr, "Wrong program usage! Try\n%s -h\n", argv[0]);
                        exit (EXIT_FAILURE);
        }
    }

    fout = (filename) ? fopen (filename, "w") : stdout;

    binary sys;
    init_system (&sys, M1, M2, eps);
    T = solar_year (sys);
    top *= T;
    dt  *= T/sqrt(2);   // for stability of the symplectic integrator

    planet omikron =
    {
        .zeta = zeta,
        .psi    = psi,
        .p_psi  = p_psi,
        .p_zeta = p_zeta
    };

    switch (option)
    {
        case 0: solver (&omikron, &sys, dt, top, fout);    break;
        case 1: solver_S8 (&omikron, &sys, dt, top, fout); break;
        case 2: adaptive_solver (&S4a, &omikron, &sys, dt, top, prec, fout); break;
        case 3: adaptive_solver (&S8, &omikron, &sys, dt, top, prec, fout); break;
        case 4: poincare (&S4a, &omikron, &sys, dt, prec, fout); break;
        default:
                fprintf (stderr, "Unknown option\n");
                exit (EXIT_FAILURE);
    }

    if (out) fclose (fout);

    return 0;
}
