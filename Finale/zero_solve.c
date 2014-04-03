#include "zero_solve.h"

double
mod (double a, double b)    // works only for positive a,b > 0
{
    double n = 0,
           r = (b == 0) ? b : b * modf (a / b, &n);

    return r;
}

double
tap (double epsilon, double phi)
{
    double E    = 1 - epsilon*epsilon,
           e    = 1 - epsilon,
           sE   = sqrt(E),
           T = 2*M_PI / (sE * sE * sE),
           t = 0;

    t = T/2 + T * M_1_PI * atan (e * tan(0.5*(phi - M_PI))/sE)
        + epsilon * sin(M_PI - phi)/(E * (1 + epsilon * cos(phi - M_PI)));

    return t;
}

double
dtap (double epsilon, double phi)
{
    double r = (1 + epsilon*cos(phi - M_PI));
    return 1.0/(r*r);
}

double
ddtap (double epsilon, double phi)
{
    double r = (1 + epsilon*cos(phi - M_PI));
    return 2.0*(epsilon)*sin(phi - M_PI)/(r * r * r);
}

void
get_phi (double * phi, double epsilon, double t)
{
    double T = 2*M_PI/pow(1 - epsilon * epsilon, 1.5),
           phi_1 = 0;

    t = mod (t, T);

    // we have to find the zero of `tap(...) - t' with Newton's method
    if (t < T/2)
    {
        do
        {
            *phi    = phi_1;
            phi_1   = *phi - (tap(epsilon, *phi) - t)/(dtap(epsilon, *phi));
        } while (fabs(phi_1 - *phi) > 1e-12 && fabs(phi_1 - *phi) != 2*M_PI);
        *phi = phi_1;
    }
    else    // function is odd over T/2, we use this to stabilize the right interval solutions
    {
        t = T - t;
        do
        {
            *phi    = phi_1;
            phi_1   = *phi - (tap(epsilon, *phi) - t)/(dtap(epsilon, *phi));
        } while (fabs(phi_1 - *phi) > 1e-12 && fabs(phi_1 - *phi) != 2*M_PI);
        *phi = 2*M_PI - phi_1;
    }
}
