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
           L      = 400,
           prec   = 1e-3;

    int option = 0,
        out    = 0,
        cont   = 0,
        arg;

    char * filename = NULL,
         * mode;

    size_t N = 100;

    FILE * fout;
    filename = NULL;
    mode = "w";
    while ((arg = getopt(argc, argv, "M:N:m:e:T:d:z:l:L:p:r:O:o:c:h")) != -1)
    {
        switch (arg)
        {
            case 'M':   M1   = atof(optarg); break;
            case 'N':   N    = atoi(optarg); break;
            case 'm':   M2   = atof(optarg); break;
            case 'e':   eps  = atof(optarg); break;
            case 'T':   top  = atof(optarg); break;
            case 'd':   dt   = atof(optarg); break;
            case 'z':   zeta = atof(optarg); break;
            case 'r':   prec = atof(optarg); break;
            case 'L':   L    = atof(optarg); break;
            case 'l':   p_psi   = atof(optarg); break;
            case 'p':   p_zeta  = atof(optarg); break;
            case 'O':   option  = atoi(optarg); break;
            case 'c':
                        filename = optarg;
                        cont     = 1;
                        mode     = "r+a";
                        break;
            case 'o':   
                        filename = optarg;
                        out      = 1;
                        break;
            case 'h':
                        fprintf (stdout, "-M:   mass of the 1st star\n");
                        fprintf (stdout, "-N:   number of points for Poincare\n");
                        fprintf (stdout, "-m:   mass of the 2nd star\n");
                        fprintf (stdout, "-e:   eccentricity of the star orbits\n");
                        fprintf (stdout, "-T:   how many star cycles we will compute\n");
                        fprintf (stdout, "-d:   time step duration\n");
                        fprintf (stdout, "-z:   zeta of the planet\n");
                        fprintf (stdout, "-l:   starting angular momentum of the planet\n");
                        fprintf (stdout, "-L:   maximum angular momentum when scanning for equilibrium points\n");
                        fprintf (stdout, "-p:   starting radial momentum of the planet\n");
                        fprintf (stdout, "-r:   precision for the adaptive step\n");
                        fprintf (stdout, "-O:   option -- various program modes\n");
                        fprintf (stdout, "   0   -> non-adaptive S4 solver\n");
                        fprintf (stdout, "   1   -> non-adaptive S8 solver\n");
                        fprintf (stdout, "   2   -> adaptive S4 solver\n");
                        fprintf (stdout, "   3   -> adaptive S8 solver\n");
                        fprintf (stdout, "   4   -> adaptive S4 solver on Poincare cross-section\n");
                        fprintf (stdout, "   5   -> adaptive S8 solver on Poincare cross-section\n");
                        fprintf (stdout, "   6   -> use bisection and adaptive S8 to find points for various L\n");
                        fprintf (stdout, "-o:   write output into a file, instead of stdout\n");
                        fprintf (stdout, "-c:   open the file to continue from where we left off\n");
                        fprintf (stdout, "-h:   print this list\n");
                        exit (EXIT_SUCCESS);
            default:
                        fprintf (stderr, "Wrong program usage! Try\n%s -h\n", argv[0]);
                        exit (EXIT_FAILURE);
        }
    }

    fprintf (stderr, "mode = %s\n", mode);
    fout = (filename) ? fopen (filename, mode) : stdout;

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
        case 3: adaptive_solver (&S8,  &omikron, &sys, dt, top, prec, fout); break;
        case 4: 
                if (cont == 0)
                    poincare (&S4a, &omikron, &sys, NULL, dt, prec, N, fout);
                else
                    Continue (&S4a, &omikron, &sys, dt, prec, N, fout);
                break;
        case 5: 
                if (cont == 0)
                    poincare (&S8, &omikron, &sys, NULL, dt, prec, N, fout);
                else
                    Continue (&S8, &omikron, &sys, dt, prec, N, fout);
                break;
        case 6:
                equilibrium (&S8, &omikron, &sys, L, dt, prec, fout);
                break;
        default:
                fprintf (stderr, "Unknown option\n");
                exit (EXIT_FAILURE);
    }

    if (out) fclose (fout);

    return 0;
}
