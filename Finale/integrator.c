#include "integrator.h"

//const double X0 =  1.35120719195966; // X0 = 1.0 / (2 - cbrt(2))
//const double X1 = -1.70241438391932; // X1 = cbrt(2) / (cbrt(2) - 1)

void
A_step (planet * omikron, 
        double dt, double c)
{
    omikron->psi    += c * dt * 0.5 * omikron->p_psi/(omikron->zeta * omikron->zeta);
    omikron->p_zeta += c * dt * 0.5 * omikron->p_psi * omikron->p_psi / (omikron->zeta * omikron->zeta * omikron->zeta);
    omikron->zeta   += c * dt * 0.5 * omikron->p_zeta;
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

    omikron->p_zeta += c * dt * 0.5*sys.rho*cos(sys.phi - omikron->psi)*dR - omikron->zeta*MR;
    omikron->p_psi  += c * dt * 0.5*sys.rho*omikron->zeta*dR*sin(sys.phi - omikron->psi);
}

void
S2 (planet * omikron,
        binary sys, double dt, double c)
{
    A_step (omikron, dt, c);
    B_step (omikron, sys, dt, c);
    A_step (omikron, dt, c);
}

void
S4 (planet * omikron,
        binary sys, double dt)
{
    S2 (omikron, sys, dt, X0);
    S2 (omikron, sys, dt, X1);
    S2 (omikron, sys, dt, X0);
}

