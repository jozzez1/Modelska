#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <getopt.h>
#include "global.h"
#include "sor.h"
#include "navier_stokes.h"

void solve (const unsigned int N,
        const unsigned int M,
        const unsigned long T,
        const unsigned long frames,
        const double Re,
        const double precision,
        double * delta)
{
    // variable declarations
    // --------------------------------------------------------
    double * tmp    = (double *) calloc (N*N, sizeof(double)),
           * zeta   = (double *) calloc (N*N, sizeof(double)),
           * swp    = (double *) calloc (N*N, sizeof(double)),
           * psi    = (double *) calloc (N*N, sizeof(double)),
           * u      = (double *) calloc (N*N, sizeof(double)),
           * v      = (double *) calloc (N*N, sizeof(double)),
           * x      = (double *) calloc (N*N, sizeof(double)),
           * y      = (double *) calloc (N*N, sizeof(double));

    point * markers = (point *) malloc (M * sizeof (point));

    char filename [15];

    double J        = cos(M_PI/N),
           norm2    = 0,
           xi       = 0,
           w        = 0;

    // solution
    // --------------------------------------------------------

    // first we set the non-zero conditions
    initial_conditions (zeta, u, N);
    init_xy (x, y, N);
    init_markers (markers, M, N);
    // we use the 1st 10 iterations to predict delta
    unsigned int i, j;
    for (i = 10; i--;)
    {
        SOR (psi, &w, &xi, &norm2, zeta, &J, precision, N);
        get_vxy_and_delta (u, v, delta, psi, N);
        swap (&zeta, &tmp, &swp);
        iterate_zeta (zeta, tmp, u, v, delta, &Re, N);
        iterate_markers (markers, M, u, v, delta, N);
        fix_zeta_boundaries (zeta, psi, N);
        printf ("i = %u\tdelta = %.2e\n", i, *delta);
    }

    // the rest we can do normally
    for (i = frames; i--;)
    {
        for (j = T; j--;)
        {
            SOR (psi, &w, &xi, &norm2, zeta, &J, precision, N);
            get_vxy (u, v, psi, N);
            swap (&zeta, &tmp, &swp);
            iterate_zeta (zeta, tmp, u, v, delta, &Re, N);
            iterate_markers (markers, M, u, v, delta, N);
            fix_zeta_boundaries (zeta, psi, N);
        }
        sprintf (filename, "anim/%06lu.jpg", frames -1 -i);
        plot_markers (markers, M, filename);
        printf ("i = %u\tdelta = %.2e\n", i, *delta);
    }
    plot_flow (x, y, u, v, N);
    system ("./anime.sh anim");

    // variable deallocation
    // --------------------------------------------------------
    if (swp)    free (swp);    // <-- could be unallocated at this point
    if (tmp)    free (tmp);    // because of that switcheroo

    free (psi);
    free (u);
    free (v);
    free (x);
    free (y);
}

int main (int argc, char ** argv)
{
    unsigned int N  = 101,
                 M  = 50,
                 arg;
    unsigned long T = 100000,
                  F = 100;
    double delta    = 1e-6,
           Re       = 1e-2,
           precision= 1e-5;

    struct option longopts[] =
    {
        {"N",       required_argument,  NULL,   'N' },
        {"M",       required_argument,  NULL,   'M' },
        {"Re",      required_argument,  NULL,   'R' },
        {"prec",    required_argument,  NULL,   'p' },
        {"time",    required_argument,  NULL,   'T' },
        {"frames",  required_argument,  NULL,   'F' },
        {"help",    no_argument,        NULL,   'h' },
        {0, 0, 0, 0}
    };

	while ((arg = getopt_long (argc, argv, "N:M:R:p:T:F:h", longopts, NULL)) != -1)
	{
        switch (arg)
        {
            case 'N': N             = atoi (optarg); break;
            case 'M': M             = atoi (optarg); break;
            case 'R': Re            = atof (optarg); break;
            case 'p': precision     = atof (optarg); break;
            case 'T': T             = atoi (optarg); break;
            case 'F': F             = atoi (optarg); break;
            case 'h':
                      printf("-N, --N:       set the matrix rank\n");
                      printf("-M, --M:       set the number of markers\n");
                      printf("-R, --Re:      Reynold's number\n");
                      printf("-p, --prec:    SOR precision condition\n");
                      printf("-T, --time:    number of time slides\n");
                      printf("-F, --frames:  number of total frames\n");
                      printf("-h, --help:    print this list\n");
                      exit (EXIT_SUCCESS);
            default:
                      printf("Wrong usage!\n");
                      printf("Try\n");
                      printf("%s -h\n", argv[0]);
                      printf("to see options.\n");
                      exit (EXIT_FAILURE);
        }
    }

    assert (N & 1);     // SOR is written in such a way, that it converges only for odd N
    solve (N, M, T, F, Re, precision, &delta);

    return 0;
}
