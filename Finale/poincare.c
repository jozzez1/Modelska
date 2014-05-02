#include "poincare.h"

unsigned int
fuzzy_compare (double a, double b)
{
    return (fabs(a - b) < 1e-6) ? 1 : 0;
}

int
poincare (scheme_fp scheme,
        planet * omikron, binary * sys, double * t_start, double dt, double precision, size_t Limit, FILE * fout)
{
    params p;
    init_params (&p);

    double t = (t_start == NULL) ? 0 : *t_start,
           z = omikron->zeta;
    size_t i = Limit;

    // set fout to be line-buffered
    setvbuf (fout, NULL, _IOLBF, 1024);
    while (i)
    {
        adaptive_step (scheme, omikron, sys, p, dt, &t, precision);
        fprintf (stderr, "z = %2.2e\n", omikron->zeta);
        if (omikron->psi <= 1e-3)
        {
            fprintf (fout, "%.12e \t%.12e \t%.12e \t%.12e \t%.12e\n",
                    t, omikron->zeta, omikron->psi, omikron->p_zeta, omikron->p_psi);
            while (omikron->psi < M_PI)
            {
                adaptive_step (scheme, omikron, sys, p, dt, &t, precision);
                if (omikron->zeta > 5*z)
                    break;
            }
            i--;
        }
        if (omikron->zeta > 100*z)
        {
            fprintf (stderr, "Stop: Planet escaped the system.\n");
            return PLANET_ESCAPE;
            break;
        }
    }
    return 0;
}

size_t
linecount (FILE * fin)
{
    size_t N = 0;
    char c;

    rewind (fin);
    do
    {
        c = fgetc (fin);
        if (c == '\n')
            N++;
    } while (c != EOF);
    rewind (fin);

    return N;
}

void
Continue (scheme_fp scheme,
        planet * omikron, binary * sys, double dt, double precision, size_t Limit, FILE * fio)
{
    size_t N = linecount (fio);
    Limit -= N;
    double t = 0;

    fprintf (stderr, "N = %lu\n", N);
    if (N)
    {
        for (size_t i = N; i--;)
            fscanf (fio, "%lf %lf %lf %lf %lf\n",
                    &t, &omikron->zeta, &omikron->psi, &omikron->p_zeta, &omikron->p_psi);
    }

    // file MUST be open to append stuff, or it will of course fail big time!
    poincare (scheme, omikron, sys, &t, dt, precision, Limit, fio);
}

double
bisect (scheme_fp scheme,
        planet * omikron, binary * sys, double Z, double dt, double precision)
{
    double z0 = omikron->zeta,  // these are two starting zeta points
           z1 = Z,              // at this zeta, p_zeta > 0
           z2 = (z0 + z1)/2,    // the midpoint
           p_z0,                // their respective radial impulses are down here
           p_z1,                // the 1st point
           p_z2;                // and the midpoint

    int signal;

    // backup of the planet settings
    planet omikron_save = { omikron->zeta, omikron->psi, omikron->p_psi, omikron->p_zeta };

    // get the p_z0
    omikron->zeta = z0;
    signal = poincare (scheme, omikron, sys, NULL, dt, precision, 2, stdout);
    if (signal == PLANET_ESCAPE) exit (EXIT_FAILURE);
    p_z0 = omikron->p_zeta;

    // restore omikron
    cpy_planets (omikron, omikron_save, 1);

    // get the p_z1
    omikron->zeta = z1;
    signal = poincare (scheme, omikron, sys, NULL, dt, precision, 2, stdout);
    if (signal == PLANET_ESCAPE) exit (EXIT_FAILURE);
    p_z1 = omikron->p_zeta;

    // again restore omikron
    cpy_planets (omikron, omikron_save, 1);

    // get p_z2
    omikron->zeta = z2;
    signal = poincare (scheme, omikron, sys, NULL, dt, precision, 2, stdout);
    if (signal == PLANET_ESCAPE) exit (EXIT_FAILURE);
    p_z2 = omikron->p_zeta;

    do
    {
        if (p_z1 * p_z0 > 0)
        {
            fprintf (stderr, "Zero is not in this interval. Stop.\n");
            exit (EXIT_FAILURE);
        }
        else if (fabs(p_z2) < 1e-8) // we found the approximate zero
        {
            z0      = z2;   // solution will be written in z0
            p_z0    = p_z2;
            break;
        }
        else if (p_z2*p_z1 > 0)
        {
            z1      = z2;
            p_z1    = p_z2;
        }
        else if (p_z2*p_z0 > 0)
        {
            z0      = z2;
            p_z0    = p_z2;
        }

        z2 = (z1 + z0) / 2;
        fprintf (stdout, "z0 = %lf\nz2 = %lf\nz1 = %lf\n",
                z0, z2, z1);
        cpy_planets (omikron, omikron_save, 1);
        omikron->zeta = z2;
        signal = poincare (scheme, omikron, sys, NULL, dt, precision, 2, stdout);
        p_z2 = omikron->p_zeta;
        if (signal == PLANET_ESCAPE) exit (EXIT_FAILURE);
    } while (fabs(z1 - z0) > 1e-6);

    return z2;
}

void
equilibrium (scheme_fp scheme,
        planet * omikron, binary * sys, double Lmax, double dt, double precision, FILE * fout)
{
    double dL   = 20,
           Z    = 40,
           z0;

    planet omikron_save =
    {
        .zeta   = omikron->zeta,
        .psi    = omikron->psi,
        .p_psi  = omikron->p_psi,
        .p_zeta = omikron->p_zeta
    };

    setvbuf (fout, NULL, _IOLBF, 1024);
    while (omikron->p_psi <= Lmax)
    {
        z0  = bisect (scheme, omikron, sys, Z, dt, precision);
        fprintf (fout, "%lf\t%lf\n",
                omikron_save.p_psi, z0);

        Z  += (z0 - omikron_save.zeta);
        omikron_save.zeta    = z0;
        omikron_save.p_psi  += dL;
        cpy_planets (omikron, omikron_save, 1);
    }
}
