#include "trotter.h"

// we pass all the needed parameters to this thingy
Trotter::Trotter (unsigned int tot, double delta_t,
        double mass_1, double mass_2,
        point3 r1, point3 r2, point3 r3,
        point3 p1, point3 p2, point3 p3)
{
    Total   = tot;
    dt      = delta_t;

    mass.push_back (mass_1);
    mass.push_back (mass_2);
    mass.push_back (1);
    
    coordinate.push_back (r1);
    coordinate.push_back (r2);
    coordinate.push_back (r3);

    momentum.push_back (p1);
    momentum.push_back (p2);
    momentum.push_back (p3);
}

const double Trotter::X0 = 1.0/(1 - cbrt(2));
const double Trotter::X1 = (-cbrt(2)) / (2 - cbrt(2));

inline void
Trotter::T_split_step (const double * c)
{
    for (unsigned int i = 3; i--;)
    {
        coordinate[i].x += *c * 0.5* dt * momentum[i].x/mass[i]; // 0.5 is from S2 scheme
        coordinate[i].y += *c * 0.5* dt * momentum[i].y/mass[i]; // so we don't have to pass it to the argument
        coordinate[i].z += *c * 0.5* dt * momentum[i].z/mass[i];
    }
}

inline void
Trotter::V_split_step (const double * c)
{
    for (unsigned int i = 3; i--;)
    {
        unsigned int ip1 = dotplus(i,1),
                     ip2 = dotplus(i,2);

        double dp1 = mass[i]*mass[ip1]/pow(pow2(coordinate[i].x - coordinate[ip1].x)
                    + pow2(coordinate[i].y - coordinate[ip1].y)
                    + pow2(coordinate[i].z - coordinate[ip1].z), 1.5),
               dp2 = mass[i]*mass[ip2]/pow(pow2(coordinate[i].x - coordinate[ip2].x)
                    + pow2(coordinate[i].y - coordinate[ip2].y)
                    + pow2(coordinate[i].z - coordinate[ip2].z), 1.5);

        momentum[i].x += *c * dt * (dp1 * (coordinate[ip1].x - coordinate[i].x) + dp2 * (coordinate[ip2].x - coordinate[i].x));
        momentum[i].y += *c * dt * (dp1 * (coordinate[ip1].y - coordinate[i].y) + dp2 * (coordinate[ip2].y - coordinate[i].y));
        momentum[i].z += *c * dt * (dp1 * (coordinate[ip1].z - coordinate[i].z) + dp2 * (coordinate[ip2].z - coordinate[i].z));
    }
}

void
Trotter::S2 (const double * c)
{
    T_split_step (c);
    V_split_step (c);
    T_split_step (c);
}

void
Trotter::S4 ()
{
    S2 (&X0);
    S2 (&X1);
    S2 (&X0);
}

void
Trotter::update_T (void)
{
    T = 0;
    for (unsigned int i = 3; i--;)
        T += (pow2(momentum[i].x) + pow2(momentum[i].y) + pow2(momentum[i].z))/mass[i];

    T /= 2;
}

void
Trotter::update_V (void)
{
    V = 0;
    V -= mass[1]/dist2(coordinate[1], coordinate[3])
        + mass[2]/dist2(coordinate[2], coordinate[3])
        + mass[1]*mass[2]/dist2(coordinate[1],coordinate[2]);
}

void
Trotter::time_step (void)
{
    S4 ();
    t += dt;
}

void
Trotter::plot_current (void)
{
    mglGraph gr (800, 800);

    mglData RX (3),
            RY (3),
            RZ (3),
            PX (3),
            PY (3),
            PZ (3);
}
