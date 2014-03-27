#include "global.h"

inline double
pow2 (const double x)
{
    return x*x;
}

inline double
mabs (const double *x)
{
    double re = *x;
    if (x < 0) re = -re;
    return re;
}

inline double
get_greater (const double *x, const double *y)
{
    double re = mabs(x),
           fe = mabs(y);
    if (re < fe) re = fe;
    return re;
}

inline void
print_array (FILE * fout, const double * xi, const unsigned long N)
{
    unsigned int i,j;
    for (i = 0; i <= N-1; i++)
    {
        for (j = 0; j <= N-1; j++)
            fprintf (fout, "% .1lf\t", xi[j + i*N]);
        fprintf (fout, "\n");
    }
    fprintf (fout, "\n");
}

void
swap (double ** A, double ** B, double ** swp)
{
    // we just swap the addresses of the arrays
    *swp = *B;
    *B = *A;
    *A = *swp;
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
