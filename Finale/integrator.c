#include "integrator.h"

void
A_step (planet * omikron, 
        double dt, double c)
{
    omikron->zeta   += c * dt * omikron->p_zeta;
}

void
B_step (planet * omikron,
        binary sys, double dt, double c)
{
    double r31 = omikron->zeta * omikron->zeta + 0.25 * sys.rho * sys.rho / (sys.M1 * sys.M1)
                    - sys.rho * omikron->zeta * cos (omikron->psi - sys.phi)/sys.M1,
           r32 = omikron->zeta * omikron->zeta + 0.25 * sys.rho * sys.rho / (sys.M2 * sys.M2)
                    + sys.rho * omikron->zeta * cos (omikron->psi - sys.phi)/sys.M2,
           dR  = 1.0/pow(r31, 1.5) - 1.0/pow(r32, 1.5),
           MR  = sys.M1/pow(r31, 1.5) + sys.M2/pow(r32, 1.5);

    double psi = omikron->psi;
    omikron->psi    += c * dt * omikron->p_psi/(omikron->zeta * omikron->zeta);
    omikron->p_zeta += c * dt * (0.5*sys.rho*cos(sys.phi - psi)*dR - omikron->zeta*MR
                        + omikron->p_psi * omikron->p_psi / (omikron->zeta * omikron->zeta * omikron->zeta));
    omikron->p_psi  += c * dt * 0.5*sys.rho*omikron->zeta*dR*sin(sys.phi - psi);

    omikron->psi = mod (omikron->psi, 2*M_PI);
}

void
S2 (planet * omikron,
        binary sys, double dt, double c)
{
    A_step (omikron, dt, c/2);
    B_step (omikron, sys, dt, c);
    A_step (omikron, dt, c/2);
}

void
S4 (planet * omikron,
        binary sys, double dt)
{
    S2 (omikron, sys, dt, X0);
    S2 (omikron, sys, dt, X1);
    S2 (omikron, sys, dt, X0);
}

void
S4a (planet * omikron,
        binary sys, params p, double dt)
{
    S4 (omikron, sys, dt);
}

void init_params (params * p)
{
    double w0 = 1 - 2*(w[1] + w[2] + w[3] + w[4] + w[5] + w[6] + w[7]);

    p->d[0] = w[7];     p->d[14] = p->d[0];
    p->d[1] = w[6];     p->d[13] = p->d[1];
    p->d[2] = w[5];     p->d[12] = p->d[2];
    p->d[3] = w[4];     p->d[11] = p->d[3];
    p->d[4] = w[3];     p->d[10] = p->d[4];
    p->d[5] = w[2];     p->d[9]  = p->d[5];
    p->d[6] = w[1];     p->d[8]  = p->d[6];
    p->d[7] = w0;

    p->c[0] = 0.5 * (w[0] + w[7]);  p->c[15] = p->c[0];
    p->c[1] = 0.5 * (w[7] + w[6]);  p->c[14] = p->c[1];
    p->c[2] = 0.5 * (w[6] + w[5]);  p->c[13] = p->c[2];
    p->c[3] = 0.5 * (w[5] + w[4]);  p->c[12] = p->c[3];
    p->c[4] = 0.5 * (w[4] + w[3]);  p->c[11] = p->c[4];
    p->c[5] = 0.5 * (w[3] + w[2]);  p->c[10] = p->c[5];
    p->c[6] = 0.5 * (w[2] + w[1]);  p->c[9]  = p->c[6];
    p->c[7] = 0.5 * (w[1] +   w0);  p->c[8]  = p->c[7];
}

void
S8 (planet * omikron,
        binary sys, params p, double dt)
{
    A_step (omikron, dt, p.c[0]);
    B_step (omikron, sys, dt, p.d[0]);
    A_step (omikron, dt, p.c[1]);
    B_step (omikron, sys, dt, p.d[1]);
    A_step (omikron, dt, p.c[2]);
    B_step (omikron, sys, dt, p.d[2]);
    A_step (omikron, dt, p.c[3]);
    B_step (omikron, sys, dt, p.d[3]);
    A_step (omikron, dt, p.c[4]);
    B_step (omikron, sys, dt, p.d[4]);
    A_step (omikron, dt, p.c[5]);
    B_step (omikron, sys, dt, p.d[5]);
    A_step (omikron, dt, p.c[6]);
    B_step (omikron, sys, dt, p.d[6]);
    A_step (omikron, dt, p.c[7]);
    B_step (omikron, sys, dt, p.d[7]);
    A_step (omikron, dt, p.c[8]);
    B_step (omikron, sys, dt, p.d[8]);
    A_step (omikron, dt, p.c[9]);
    B_step (omikron, sys, dt, p.d[9]);
    A_step (omikron, dt, p.c[10]);
    B_step (omikron, sys, dt, p.d[10]);
    A_step (omikron, dt, p.c[11]);
    B_step (omikron, sys, dt, p.d[11]);
    A_step (omikron, dt, p.c[12]);
    B_step (omikron, sys, dt, p.d[12]);
    A_step (omikron, dt, p.c[13]);
    B_step (omikron, sys, dt, p.d[13]);
    A_step (omikron, dt, p.c[14]);
    B_step (omikron, sys, dt, p.d[14]);
    A_step (omikron, dt, p.c[15]);
}

void
Poisson (planet * deriv,
        planet * omikron, binary sys)
{
    double r31 = omikron->zeta * omikron->zeta + 0.25 * sys.rho * sys.rho / (sys.M1 * sys.M1)
                    - sys.rho * omikron->zeta * cos (omikron->psi - sys.phi)/sys.M1,
           r32 = omikron->zeta * omikron->zeta + 0.25 * sys.rho * sys.rho / (sys.M2 * sys.M2)
                    + sys.rho * omikron->zeta * cos (omikron->psi - sys.phi)/sys.M2,
           dR  = 1.0/pow(r31, 1.5) - 1.0/pow(r32, 1.5),
           MR  = sys.M1/pow(r31, 1.5) + sys.M2/pow(r32, 1.5);

    deriv->zeta     = omikron->p_zeta;
    deriv->psi      = omikron->p_psi / (omikron->zeta * omikron->zeta);
    deriv->p_zeta   = omikron->p_psi * omikron->p_psi / (omikron->zeta * omikron->zeta * omikron->zeta)
                        + 0.5*sys.rho * cos(sys.phi - omikron->psi)*dR - omikron->zeta*MR;
    deriv->p_psi    = 0.5*sys.rho * omikron->zeta * dR * sin(sys.phi - omikron->psi);
}

void
sum_planets (planet * omikron, planet pluto, double c)
{
    omikron->zeta   += c * pluto.zeta;
    omikron->psi    += c * pluto.psi;
    omikron->p_zeta += c * pluto.p_zeta;
    omikron->p_psi  += c * pluto.p_psi;
}

void
RK4 (planet * omikron,
        binary * sys, double t, double dt)
{
    planet k1, k2, k3, k4, y;
    sum_planets (&y, *omikron, 1);

    get_position (sys, t);
    Poisson (&k1, &y, *sys);

    get_position (sys, t + 0.5*dt);
    sum_planets (&y, k1, dt*0.5);
    Poisson (&k2, &y, *sys);

    sum_planets (&y, *omikron, 1);
    sum_planets (&y, k2, dt*0.5);
    Poisson (&k3, &y, *sys);

    sum_planets (&y, *omikron, 1);
    get_position (sys, t + dt);
    sum_planets (&y, k3, dt);
    Poisson (&k4, &y, *sys);

    sum_planets (omikron, k1, dt*1.0/6);
    sum_planets (omikron, k2, dt*1.0/3);
    sum_planets (omikron, k3, dt*1.0/3);
    sum_planets (omikron, k4, dt*1.0/6);
}

void
adaptive_step (void (* scheme) (planet*, binary, params, double),
        planet * omikron, binary * sys, params p, double dt,
        double * t, double precision)
{
    double t_old = *t,
           error = 0;

    planet omikron_old = *omikron,
           omikron_prev= *omikron;

    // we have to make the 1st step here before we do anything else
    get_position (sys, *t);
    scheme (omikron, *sys, p, dt);

    unsigned int i = 1,
                 k = 2;
    do
    {
        *t = t_old;
        omikron_prev = *omikron;
        *omikron = omikron_old;

        for (i = k; i--;)
        {
            get_position (sys, *t);
            scheme (omikron, *sys, p, dt/k);
            *t += dt/k;
        }

        error = (omikron->zeta - omikron_prev.zeta) * (omikron->zeta - omikron_prev.zeta)
              + (omikron->psi - omikron_prev.psi) * (omikron->psi - omikron_prev.psi)
              + (omikron->p_psi - omikron_prev.p_psi) * (omikron->p_psi - omikron_prev.p_psi)
              + (omikron->p_zeta - omikron_prev.p_zeta) * (omikron->p_zeta - omikron_prev.p_zeta);
        error = sqrt (error);
//        k++;
        k *= 2; // we have to use Romberg sequence, or symplectic integrator will fail
    } while (fabs(error) > precision);
    *t = t_old + dt;
    get_position (sys, t_old);
}

void
solver (planet * omikron, binary * sys,
        double dt, double T, FILE * fout)
{
    double t = 0;   // of course
    while (t < T)
    {
        get_position (sys, t);
        S4 (omikron, *sys, dt);
        t += dt;
        fprintf (fout, "%.12lf \t%.12lf \t%.12lf \t%.12lf \t%.18lf\n",
                t, sys->rho, sys->phi, omikron->zeta, omikron->psi);
    }
}

void
solver_S8 (planet * omikron, binary * sys,
        double dt, double T, FILE * fout)
{
    params p;
    init_params (&p);

    double t = 0;
    while (t < T)
    {
        get_position (sys, t);
        S8 (omikron, *sys, p, dt);
        t += dt;
        fprintf (fout, "%.12lf \t%.12lf \t%.12lf \t%.12lf \t%.18lf\n",
                t, sys->rho, sys->phi, omikron->zeta, omikron->psi);
    }
}

void
solver_RK4 (planet * omikron, binary * sys,
        double dt, double T, FILE * fout)
{
    double t = 0;
    while (t < T)
    {
        RK4 (omikron, sys, t, dt);
        t += dt;
        fprintf (fout, "%.12lf \t%.12lf \t%.12lf \t%.12lf \t%.18lf\n",
                t, sys->rho, sys->phi, omikron->zeta, omikron->psi);
    }
}

void
adaptive_solver (void (*scheme) (planet *, binary, params, double),
        planet * omikron, binary * sys, double dt,
        double T, double precision, FILE * fout)
{
    params p;
    init_params (&p);

    double t = 0;
    while (t < T)
    {
        adaptive_step (scheme, omikron, sys, p, dt, &t, precision);
        fprintf (fout, "%.12lf \t%.12lf \t%.12lf \t%.12lf \t%.18lf\n",
                t, sys->rho, sys->phi, omikron->zeta, omikron->psi);
    }
}
