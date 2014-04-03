#include "star_system.h"

void
init_system (binary * u, double M1, double M2, double p_phi, double p_rho)
{
    double mu       = M1 * M2 /(M1 + M2),
           alpha    = 8 * mu * mu * M1 * M2;

    // now we check if p_phi and p_rho are valid
    double Limit_p_phi  = sqrt(4*M1 * M2),
           limit_p_phi  = sqrt(4*M1 * M2 * (1 - sqrt((M1*M1 + M2*M2)/((M1 + M2)*(M1 + M2)))));
    
    assert (p_phi < Limit_p_phi);
    assert (p_phi >= limit_p_phi);

    double p_phi2       = p_phi * p_phi,
           Limit_p_rho  = sqrt(4*M1*M2 - p_phi2),
           limit_p_rho  = sqrt(4*M1*M2 - p_phi2 - alpha/p_phi2);

    assert (p_rho < Limit_p_rho);
    assert (p_rho >= limit_p_rho);

    // so far so good, all the parameters were accepted
    u->M1       = M1;
    u->M2       = M2;
    u->mu       = mu;
    u->alpha    = alpha;
    u->p_phi    = p_phi;
    u->p_rho    = p_rho;
    u->E        = 0.5*(p_phi2 + p_rho*p_rho) - 2*M1*M2;
    u->epsilon  = sqrt(1 + 2*u->E*p_phi2/(alpha*alpha));
}

void
trajectory (double * rho, double * phi,
        binary sys, double t)
{
    t *= (sys.alpha * sys.alpha) / (sys.p_phi * sys.p_phi * sys.p_phi);
    get_phi (phi, sys.epsilon, t);
    *rho = (sys.p_phi * sys.p_phi) / (sys.alpha * (1 + sys.epsilon*cos(*phi)));
}


