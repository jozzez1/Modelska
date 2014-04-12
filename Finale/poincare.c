#include "poincare.h"

void
poincare (void (*scheme) (planet *, binary, params, double),
        planet * omikron, binary * sys, double dt, double precision, size_t Limit, FILE * fout)
{
    params p;
    init_params (&p);

    double t = 0,
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
        if (omikron->zeta > 5*z)
        {
            fprintf (stderr, "Stop: Planet escaped the system.\n");
            break;
        }
    }
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
Continue (void (*scheme) (planet *, binary, params, double),
        planet * omikron, binary * sys, double dt, double precision, size_t Limit, FILE * fio)
{
    size_t N = linecount (fio);
    double t = 0;

    fprintf (stderr, "N = %lu\n", N);
    if (N)
    {
        for (size_t i = N; i--;)
        {
            fscanf (fio, "%lf %lf %lf %lf %lf\n",
                    &t, &omikron->zeta, &omikron->psi, &omikron->p_zeta, &omikron->p_psi);
        }
    }

    // file MUST be open to append stuff, or it will of course fail big time!
    poincare (scheme, omikron, sys, dt, precision, Limit, fio);
}

