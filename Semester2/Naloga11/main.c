#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <mgl2/data.h>
#include <mgl2/mgl.h>

static inline double
pow2 (double x) { return x*x; };

void
tridiag_params (gsl_vector * main_diag, gsl_vector * values,
        gsl_vector * f1, gsl_vector * f0, double ht)
{
    unsigned int N = f1->size,
                 i;

    double f10          = gsl_vector_get (f1,0),
           differenceV  = gsl_vector_get(f1,1) - gsl_vector_get (f0, 1),
           differenceM  = gsl_vector_get (f1, 0) - gsl_vector_get (f1, 2),
           h            = 1.0 / N,
           hht2         = pow2 (h / ht);
    gsl_vector_set (values, 0, (-1)*hht2 * pow2(differenceV) - h*sin (f10));
    gsl_vector_set (main_diag, 0, -1 -0.25*pow2(differenceM));

    for (i = N-3; i--;)
    {
        differenceV = gsl_vector_get (f1, i+2) - gsl_vector_get (f0, i+2);
        differenceM = gsl_vector_get (f1, i+1) - gsl_vector_get (f1, i+3);
        gsl_vector_set (values, i+1, (-1)*hht2 * pow2(differenceV));
        gsl_vector_set (main_diag, i+1, -2 -0.25*pow2(differenceM));
    }
}

void
Fnew (gsl_vector * solution, gsl_vector * main_diag, gsl_vector * values,
        gsl_vector * f1, gsl_vector * f0, gsl_vector * side_diag, double ht)
{
    unsigned int N = f0->size;
    tridiag_params (main_diag, values, f1, f0, ht);
    gsl_vector_view subsolution = gsl_vector_subvector (solution, 1, N-2);

    // now we solve this system
    gsl_linalg_solve_symm_tridiag (main_diag, side_diag, values, &subsolution.vector);
    gsl_vector_set (solution, N-1, 0);
    gsl_vector_set (solution, 0, sin(gsl_vector_get (f1, 0))/N - gsl_vector_get (solution, 1));
}

void
fnew (gsl_vector * f_updated, gsl_vector * F, gsl_vector * f1, gsl_vector * f0, double ht)
{
    unsigned int N = f0->size, i;
    double hth2 = pow2 (ht*N);

    for (i = N-2; i--;)
    {
        double f11 = gsl_vector_get (f1, i+2),
               f1n = gsl_vector_get (f1, i+1),
               f10 = gsl_vector_get (f1, i),
               f00 = gsl_vector_get (f0, i+1),
               F1  = gsl_vector_get (F, i+2),
               F0  = gsl_vector_get (F, i),
               Fn  = gsl_vector_get (F, i+1),
               com = hth2 * (0.5*(f11 - f10)*(F1 - F0) + Fn*(f11 + f10 - 2*f1n)) + 2*f1n - f00;

        gsl_vector_set (f_updated, i+1, com);
    }
    double fT2 = gsl_vector_get (f_updated, 1),
           FT1 = gsl_vector_get (F, 1),
           fT3 = gsl_vector_get (f_updated, 2),
           fTn = gsl_vector_get (f_updated, N-3),
           fTN = 2*gsl_vector_get (f_updated, N-2) - fTn,
           minus_fT0 = 2*cos(fT2)/(N * FT1) + fT3;

    gsl_vector_set (f_updated, 0, minus_fT0);
    gsl_vector_set (f_updated, N-1, fTN);
}

void
time_step (gsl_vector * FT, gsl_vector * fT, gsl_vector * main_diag, gsl_vector * values,
        gsl_vector * F, gsl_vector * side_diag,
        gsl_vector * f1, gsl_vector * f0, double ht)
{
    fnew (fT, F, f1, f0, ht);
    Fnew (FT, main_diag, values, fT, f1, side_diag, ht);
}

void
get_xy (HMDT x, HMDT y, gsl_vector * fT, unsigned int N)
{
    mgl_data_set_value (x, 0, 0, 0, 0);
    mgl_data_set_value (y, 0, 0, 0, 0);

    unsigned int i;
    for (i = N-1; i--;)
    {
        mgl_data_set_value (x, cos(gsl_vector_get (fT, N-2-i))/N + mgl_data_get_value (x, N-2-i, 0, 0), N-1-i, 0, 0);
        mgl_data_set_value (y, sin(gsl_vector_get (fT, N-2-i))/N + mgl_data_get_value (y, N-2-i, 0, 0), N-1-i, 0, 0);
    }
    // these vectors have one extra component
    mgl_data_set_value (x, cos(gsl_vector_get (fT, N-1))/N + mgl_data_get_value (x, N-1, 0, 0), N, 0, 0);
    mgl_data_set_value (y, sin(gsl_vector_get (fT, N-1))/N + mgl_data_get_value (y, N-1, 0, 0), N, 0, 0);
}

void
get_vxy (HMDT X, HMDT Y, HMDT VX, HMDT VY, gsl_vector * vx, gsl_vector * vy,
        HMDT x0, HMDT y0, HMDT x1, HMDT y1, int N, double ht)
{
    unsigned int i;
    for (i = N+1; i--;)
    {
        gsl_vector_set (vx, i, (mgl_data_get_value (x1, i, 0, 0) - mgl_data_get_value (x0, i, 0, 0))/ht);
        gsl_vector_set (vy, i, (mgl_data_get_value (y1, i, 0, 0) - mgl_data_get_value (y0, i, 0, 0))/ht);
    }

    // and now we set the rest
    for (i = 10; i--;)
    {
        mgl_data_set_value (X, mgl_data_get_value (x1, (i+1)*N/10, 0, 0), i, 0, 0);
        mgl_data_set_value (Y, mgl_data_get_value (y1, (i+1)*N/10, 0, 0), i, 0, 0);
        mgl_data_set_value (VX, gsl_vector_get (vx, (i+1)*N/10), i, 0, 0);
        mgl_data_set_value (VY, gsl_vector_get (vy, (i+1)*N/10), i, 0, 0);
    }
}

void
plot_xy (HMDT x, HMDT y, gsl_vector * FT,
        HMDT X, HMDT Y, HMDT VX, HMDT VY, char * filename)
{
    // we use FT for colors
    HMDT colors = mgl_create_data ();
    mgl_data_set_vector (colors, FT);

    HMGL gr = mgl_create_graph (800, 400);
    mgl_set_range_val (gr, 'x', 1.2, -1.2);     // we invert the 'x' axis
    mgl_set_range_val (gr, 'y', 1.2, 0);         // we invert the 'y' axis
    mgl_set_origin (gr, 1.2, 0, 0);
    mgl_axis (gr, "xy", "", "");
    mgl_tens_xy (gr, x, y, colors, "", "");
    mgl_traj_xy (gr, X, Y, VX, VY, "", "");
    mgl_write_frame (gr, filename, "");
    mgl_delete_graph (gr);
}

void
solve (unsigned int N, unsigned long int T, unsigned int frames, double ht, double fi0)
{
    // some variable declarations
    // -------------------------------------------------
    gsl_vector * f0         = gsl_vector_alloc (N),
               * f1         = gsl_vector_alloc (N),
               * FT         = gsl_vector_alloc (N),
               * fT         = gsl_vector_alloc (N),
               * F          = gsl_vector_alloc (N),
               * main_diag  = gsl_vector_alloc (N-2),
               * side_diag  = gsl_vector_alloc (N-3),
               * values     = gsl_vector_alloc (N-2),
               * vx         = gsl_vector_alloc (N+1),
               * vy         = gsl_vector_alloc (N+1);

    HMDT x0 = mgl_create_data_size (N+1, 0, 0),
         y0 = mgl_create_data_size (N+1, 0, 0),
         x1 = mgl_create_data_size (N+1, 0, 0),
         y1 = mgl_create_data_size (N+1, 0, 0),
         VX = mgl_create_data_size (10, 0, 0),
         VY = mgl_create_data_size (10, 0, 0),
         X  = mgl_create_data_size (10, 0, 0),
         Y  = mgl_create_data_size (10, 0, 0);

    char filename[15];

    // here we start solving
    // -------------------------------------------------
    gsl_vector_set_all (f0, fi0);
    gsl_vector_set_all (f1, fi0);
    gsl_vector_set_all (side_diag, 1);
    Fnew (FT, main_diag, values, f1, f0, side_diag, ht);

    unsigned long int i,j;
    for (i = frames; i--; )
    {
        // we will skip T frames
        for (j = 2*T; j--; )
        {
            gsl_vector_memcpy (F, FT);
            time_step (FT, fT, main_diag, values, F, side_diag, f1, f0, ht);
            gsl_vector_memcpy (f0, f1);
            gsl_vector_memcpy (f1, fT);
        }

        // and only plot then
        sprintf (filename, "anim/%06lu.jpg", frames-1-i);
        get_xy (x0, y0, f0, N);
        get_xy (x1, y1, f1, N);
        get_vxy (X, Y, VX, VY, vx, vy, x0, y0, x1, y1, N, ht);
        plot_xy (x1, y1, FT, X, Y, VX, VY, filename);
    }

    // free all those variables
    // -------------------------------------------------
    gsl_vector_free (f0);
    gsl_vector_free (f1);
    gsl_vector_free (FT);
    gsl_vector_free (fT);
    gsl_vector_free (F);
    gsl_vector_free (side_diag);
    gsl_vector_free (main_diag);
    gsl_vector_free (values);
    gsl_vector_free (vx);
    gsl_vector_free (vy);
    mgl_delete_data (x0);
    mgl_delete_data (y0);
    mgl_delete_data (x1);
    mgl_delete_data (y1);
    mgl_delete_data (VX);
    mgl_delete_data (VY);
    mgl_delete_data (X);
    mgl_delete_data (Y);

    // use mencoder to make animation from those pics
    // -------------------------------------------------
    system ("./anime.sh anim");
}

int main (int argc, char ** argv)
{
    unsigned int N      = 500,
                 frames = 2000;
    unsigned long int T = 10000;    // number of frames we skip
    double ht           = 1e-5,
           fi0          = M_PI*0.5;

    if (argc != 2)
    {
        printf ("Wrong program usage!\n");
        printf ("Try:\n");
        printf ("%s <fi0 [units of pi]>\n", argv[0]);
        exit (EXIT_FAILURE);
    }

    fi0 += M_PI * atof (argv[1]);

    solve (N, T, frames, ht, fi0);
    exit (EXIT_SUCCESS);
}
