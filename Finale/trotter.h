#ifndef TROTTER_H
#define TROTTER_H

#include <iostream>
#include <vector>
#include <mgl2/data.h>
#include <mgl2/mgl.h>
#include "global.h"

typedef std::vector <point3> matrix;
typedef std::vector <double> vector;

class Trotter
{
    public:
        // constructor
        Trotter (unsigned int, double,
                double, double,
                point3, point3, point3,
                point3, point3, point3);

        // read arguments from a file
        Trotter (std::ifstream);

        // destructor
        ~Trotter (void) { };

        double energy (void)    { return T + V; };
        double potential (void) { return V; };
        double kinetic (void)   { return T; };

        void time_step (void);
        void plot_current (void);
        void plot_orbits (void);

    private:
        matrix coordinate;
        matrix momentum;
        vector mass;

        double T;   // kinetic energy
        double V;   // potential energy
        double dt;  // time step "length"
        double t;   // current time

        unsigned int Total; // number of iterations

        // we have two constants for S4
        static const double X0;
        static const double X1;

        // propagator "components"
        void V_split_step (const double *);
        void T_split_step (const double *);

        // full split-step schemes
        void S2 (const double *);
        void S4 (void);

        // energy and stuff
        void update_T (void);
        void update_V (void);
};

#endif
