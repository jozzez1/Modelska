#include "star_system.h"

void
init_system (binary * u, double M1, double M2, double epsilon)
{
    double mu       = M1 * M2 /(M1 + M2),
           alpha    = 8 * mu * mu * M1 * M2;

    // so far so good, all the parameters were accepted
    u->M1       = M1;
    u->M2       = M2;
    u->mu       = mu;
    u->alpha    = alpha;
    u->p_phi    = sqrt(2*mu*alpha*(1 - epsilon)); //p_phi;
    u->epsilon  = epsilon;
}

void
trajectory (double * rho, double * phi,
        binary sys, double t)
{
    t *= (sys.alpha * sys.alpha) / (sys.p_phi * sys.p_phi * sys.p_phi);
    get_phi (phi, sys.epsilon, t);
    *rho = (2*sys.mu * (1 - sys.epsilon) / (1 + sys.epsilon*cos(*phi - M_PI)));
}

double
solar_year (binary sys)
{
    double T = 2*M_PI / pow(1 - sys.epsilon * sys.epsilon, 1.5);
    return T * sys.p_phi * sys.p_phi * sys.p_phi/(sys.alpha * sys.alpha);
}

