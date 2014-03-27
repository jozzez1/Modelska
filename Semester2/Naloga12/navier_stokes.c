#include "navier_stokes.h"

void
initial_conditions (double * zeta, double * u,
        const unsigned int N)
{
    unsigned int i;
    for (i = N; i--;)
    {
        // the lowest edge for u is a bit different
        u [i + (N-1)*N] = 1;

        // zeta(i,j) one before last has to be fixed
        zeta [i + (N-2)*N] = u[i + (N-1)*N] * (double) N;

        // and the very last line -- the velocity for psi
        zeta [i + (N-1)*N] = (-2.0 * N);
    }
} 

void
get_vxy_and_delta (double * u, double * v, double * delta,
        const double * psi, const unsigned int N)
{
    unsigned int i,j,k;
    double over_two_h = 0.5 * N,
           vxmax      = 0,
           vymax      = 0;
    for (i = N-2; i--;)
    {
        for (j = N-2; j--;)
        {
            k = j+1 + (i+1)*N;
            u [k] = (psi[k + N] - psi[k - N]) * over_two_h; 
            v [k] = (psi[k - 1] - psi[k + 1]) * over_two_h;

            vxmax = get_greater (&u[k], &vxmax);
            vymax = get_greater (&v[k], &vymax);
        }
    }
    // and now we correct the time step -- it's 10-times smaller than the limit one
    *delta = 1.0 / (40 * N * sqrt(pow2(vxmax) + pow2(vymax)));
}

void
get_vxy (double * u, double * v,
        const double * psi, const unsigned int N)
{
    unsigned int i,j,k;
    double over_two_h = 0.5 * N;
    for (i = N-2; i--;)
    {
        for (j = N-2; j--;)
        {
            k = j+1 + (i+1)*N;
            u [k] = (psi[k - N] - psi[k + N]) * over_two_h; // we have to change the sign
            v [k] = (psi[k + 1] - psi[k - 1]) * over_two_h; // because our coordinate system is inverted
        }
    }
}

void
iterate_zeta (double * zeta,
        const double * tmp_zeta, const double * u, const double * v, const double * delta,  const double * Re, const unsigned int N)
{
    unsigned int i,j,k;
    double a = *delta * (double) N,
           b = ((double) N) / *Re;
    for (i = N-2; i--;)
    {
        for (j = N-2; j--;)
        {
            k = j+1 + (i+1)*N;
            zeta[k] = tmp_zeta[k] + a*(b*(tmp_zeta[k+1] + tmp_zeta[k-1] + tmp_zeta[k+N] + tmp_zeta[k-N] - 4*tmp_zeta[k])
                    - 0.5*(u[k+1]*tmp_zeta[k+1] - u[k-1]*tmp_zeta[k-1] + v[k+N]*tmp_zeta[k+N] - v[k-N]*tmp_zeta[k-N]));
        }
    }
}

void
fix_zeta_boundaries (double * zeta, const double * psi, const unsigned int N)
{
    unsigned int k;
    double over_h = (double) N;
    for (k = N; k--;)
    {
        // the upper boundary: zeta_ij, i = 0, j = 0,N-1
        zeta [k] = 2*pow2(N) * psi[k + N];

        // left boundary
        zeta [k*N] = 2*pow2(N) * psi[1 + k*N];

        // right boundary
        zeta [N-1 + k*N] = 2*pow2(N) * psi[N-2 + k*N];

        // bottom boundary, yay! velocity is non-zero
        zeta [k + (N-1)*N] = 2*over_h * (over_h * psi[k + (N-2)*N] - 1);
    }
}

void
init_xy (double * x, double * y, const unsigned int N)
{
    unsigned int i,j,k;
    for (i = N; i--;)
    {
        for (j = N; j--;)
        {
            k = i + j*N;    // MathGL requires its arguments to be in col. major format
            x[k] = i*1.0/(N-1);
            y[k] = j*1.0/(N-1);
            printf ("(x,y) = (%.1lf,%.1lf)\n", x[k], y[k]);
        }
    }
}

void
plot_flow (const double * x, const double * y, const double * u, const double * v, const unsigned int N)
{
    HMDT X = mgl_create_data_size (N, N, 1),
         Y = mgl_create_data_size (N, N, 1),
         VX= mgl_create_data_size (N, N, 1),
         VY= mgl_create_data_size (N, N, 1);

    mgl_data_set_double (X, x, N, N, 1);
    mgl_data_set_double (Y, y, N, N, 1);
    mgl_data_set_double (VX,u, N, N, 1);
    mgl_data_set_double (VY,v, N, N, 1);

    printf ("X=%lf x=%lf\n", mgl_data_get_value (X, 4, 2, 0), x[4 + 2*N]);

    HMGL gr = mgl_create_graph (800, 800);
    mgl_set_quality (gr, 6);
    mgl_set_range_val (gr, 'x', 0, 1);
    mgl_set_range_val (gr, 'y', 1, 0);
    mgl_set_origin (gr, 0, 0, 0);
    mgl_axis (gr, "xy", "", "");
    mgl_box (gr);
    mgl_flow_xy (gr, X, Y, VX,VY, "v", "");
    mgl_write_png (gr, "vect.png", "Velocity field");

    mgl_delete_graph (gr);
    mgl_delete_data (VX);
    mgl_delete_data (VY);
    mgl_delete_data (X);
    mgl_delete_data (Y);
}

void
init_markers (point * markers, const unsigned M, const unsigned N)
{
    double h = 1 /(2 + sqrt(M)),
           x = 1-h,
           y = 1-h;
    unsigned int i = M-1;
    while (y >= h && i)
    {
        x = 1-h;
        while (x >= h && i)
        {
            markers[i].x = x;
            markers[i].y = y;

            i--;
            x -= h;
        }
        y -= h;
    }
    assert (!i);
}

void
iterate_markers (point * markers,
        const unsigned int M, const double * u, const double * v, const double * delta, const unsigned int N)
{
    unsigned int i;
    for (i= M; i--;)
    {
        unsigned int k = (unsigned int) N * markers[i].y + N*N * (unsigned int) markers[i].x;
        markers[i].x += (*delta)*u[k];
        markers[i].y += (*delta)*v[k];
    }
}

void
plot_markers (point * markers, const unsigned int M, const char * filename)
{
    HMDT x = mgl_create_data_size (M, 0, 0),
         y = mgl_create_data_size (M, 0, 0);

    unsigned int i;
    for (i = M; i--;)
    {
        mgl_data_set_value (x, markers[i].x, i, 0, 0);
        mgl_data_set_value (y, markers[i].y, i, 0, 0);
    }

    HMGL gr = mgl_create_graph (800, 800);
    mgl_set_quality (gr, 6);
    mgl_set_range_val (gr, 'x', 0, 1);
    mgl_set_range_val (gr, 'y', 1, 0);
    mgl_set_origin (gr, 0, 0, 0);
    mgl_axis (gr, "xy", "", "");
    mgl_box (gr);
    mgl_plot_xy (gr, x, y, " #.", "");
    mgl_write_png (gr, filename, "Position of the markers");

    mgl_delete_graph (gr);
    mgl_delete_data (x);
    mgl_delete_data (y);
}
