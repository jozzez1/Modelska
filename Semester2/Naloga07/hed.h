// hed.h
///////////
#ifndef __HEADER_MOD207
#define __HEADER_MOD207

#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_linalg.h>

typedef struct
{
	double x,	// x-coordinate of the vertex
	       y,	// y-coordinate of the vertex
	       val;	// velocity at that point

	int attribute;	// for the program `Triangle' (edge or not)
} zax;

typedef struct
{
	size_t T,	// number of triangles
	       N;	// number of vertices

	zax * v;	// velocity at each vertex

	double ** No;
	int ** To;
	
	double *** A,
	       ** skal,
	       ** S,
	       * g,
	       * c;
} tr;

void destroy_tr (tr * u)
{
	return;
	int i, j;
	// free velocities and vectors
	if (u->v) free (u->v);
	if (u->g) free (u->g);
	if (u->c) free (u->c);

	// now the matrices
	for (i = 0; i <= u->N-1; i++)
	{
		if (u->No[i]) free (u->No[i]);
		if (u->S[i]) free (u->S[i]);
	}
	if (u->No) free (u->No);
	if (u->S) free (u->S);

	for (i = 0; i <= u->T-1; i++)
	{
		if (u->To[i]) free (u->To[i]);
		if (u->skal[i]) free (u->skal[i]);

		for (j = 0; j <= 2; j++)
			if (u->A[i][j]) free (u->A[i][j]);

		if (u->A[i]) free (u->A[i]);
	}
	if (u->To) free (u->To);
	if (u->skal) free (u->skal);
	if (u->A) free (u->A);

	// and lastly free the pointer
	if (u) free (u);
}

// creates the "*.poly" file
void PolyFile (int segments)
{
	FILE * fout;
	int N, i;
	if (segments)
	{
		fout = fopen ("semicircle.poly", "w");
		N = segments;
	}
	else
	{
		fout = fopen ("batman.poly", "w");
		N = 6;
	}

	// create edges
	zax * edge = (zax *) malloc (N * sizeof (zax));

	if (!segments)
	{
		edge [0].x =  0; edge[0].y =  0; edge[0].val = 0; edge[0].attribute = 2;
		edge [1].x = 15; edge[1].y =  0; edge[1].val = 1; edge[1].attribute = 2;
		edge [2].x = 15; edge[2].y = 15; edge[2].val = 2; edge[2].attribute = 2;
		edge [3].x = 10; edge[3].y = 10; edge[3].val = 3; edge[3].attribute = 2;
		edge [4].x =  5; edge[4].y = 10; edge[4].val = 4; edge[4].attribute = 2;
		edge [5].x =  0; edge[5].y = 15; edge[5].val = 5; edge[5].attribute = 2;
	}
	else
	{
		// radius = 1
		for (i = 0; i <= N-1; i++)
		{
			edge[i].x	= cos(M_PI*i/(N-1));
			edge[i].y	= sin(M_PI*i/(N-1));
			edge[i].val	= i;
			edge[i].attribute= 2;
		}
	}

	fprintf (fout, "# File automatically created\n");
	fprintf (fout, "%d %d %d %d\n",
			N, 2, 0, 1);
	
	fprintf (fout, "# Boundary vertices\n");
	for (i = 0; i <= N-1; i++)
		fprintf (fout, "%d\t %e\t %e\t %d\n",
				(int) edge[i].val, edge[i].x, edge[i].y, edge[i].attribute);
	
	fprintf (fout, "# Segments\n");
	fprintf (fout, "%d 1\n", N);
	for (i = 0; i <= N-1; i++)
		fprintf (fout, "%d\t %d\t %d\t %d\n",
				i, i, (i+1)%N, edge[i].attribute);

	fprintf (fout, "# No holes\n");
	fprintf (fout, "0\n");

	fclose (fout);
	free (edge);
}

void getTriangles (int segments, tr * u)
{
	int N, i;
	FILE * fnode,
	     * fele;

	if (segments != 0)
	{
		fnode = fopen ("semicircle.1.node", "r");
		fele  = fopen ("semicircle.1.ele", "r");
	}
	else
	{
		fnode = fopen ("batman.1.node", "r");
		fele  = fopen ("batman.1.ele", "r");
	}

	if (!fnode || !fele)
	{
		printf ("No input. Exiting.\n");
		exit (EXIT_FAILURE);
	}

	// let's do the nodes first
	fscanf (fnode, "%d %*d %*d %*d\n", &N);

	u->N = N;
	u->No = (double **) malloc (N * sizeof (double *));
	for (i = 0; i <= N-1; i++)
	{
		u->No[i] = (double *) malloc (4 * sizeof (double));
		fscanf (fnode, "%lf\t%lf\t%lf\t%lf\n",
				&u->No[i][0],	// vertice index
				&u->No[i][1],	// "x" component
				&u->No[i][2],	// "y" component
				&u->No[i][3]);	// if we are on the edge or not
	}
	fclose (fnode);

	// an now the triangles
	fscanf (fele, "%d %*d %*d\n", &N);

	u->T = N;
	u->To = (int **) malloc (N * sizeof (int *));
	for (i = 0; i <= N-1; i++)
	{
		u->To[i] = (int *) malloc (4 * sizeof (int));
		fscanf (fele, "%d %d %d %d\n",
				&u->To[i][0],	// triangle index
				&u->To[i][1],	// 1st vertice index
				&u->To[i][2],	// 2nd vertice index
				&u->To[i][3]);	// 3rd vertice index
	}
	fclose (fele);
}

// returns surface of `i-th' the triangle
double surf (tr * u, int i)
{
	double x [3],
	       y [3],
	       T;

	x[0] = u->No[u->To[i][1]][1];
	x[1] = u->No[u->To[i][2]][1];
	x[2] = u->No[u->To[i][3]][1];

	y[0] = u->No[u->To[i][1]][2];
	y[1] = u->No[u->To[i][2]][2];
	y[2] = u->No[u->To[i][3]][2];

	T  = 0;
	T += x[0]*y[1] + y[0]*x[2] + x[1]*y[2];
	T -= y[0]*x[1] + x[0]*y[2] + y[1]*x[2];
	T /= 2;
	T = fabs (T);

	return T;
}

// double we gotta get all the matrices 'A'
void calculateAscal (tr * u)
{
	int t, m, n;

	u->A = (double ***) malloc (u->T * sizeof (double **));
	u->skal = (double **) malloc (u->T * sizeof (double *));
	for (t = 0; t <= u->T-1; t++)
	{
		u->skal[t] = (double *) malloc (3 * sizeof (double));
		u->A[t] = (double **) malloc (3 * sizeof (double *));
		for (m = 0; m <= 2; m++)
			u->A[t][m] = (double *) malloc (3 * sizeof (double));
	}

	for (t = 0; t <= u->T-1; t++)
	{
		double T = surf (u, t);
		for (m = 0; m <= 2; m++)
		{
			double xm [3],
			       ym [3],
			       wm = u->No[u->To[t][(m+0)%3 + 1]][3];

			if (wm == 0)
				wm = 1;
			else if (wm > 0)
				wm = 0;

			xm[0] = u->No[u->To[t][(m+0)%3 + 1]][1];
			xm[1] = u->No[u->To[t][(m+1)%3 + 1]][1];
			xm[2] = u->No[u->To[t][(m+2)%3 + 1]][1];

			ym[0] = u->No[u->To[t][(m+0)%3 + 1]][2];
			ym[1] = u->No[u->To[t][(m+1)%3 + 1]][2];
			ym[2] = u->No[u->To[t][(m+2)%3 + 1]][2];

			for (n = 0; n <= 2; n++)
			{
				double xn [3],
				       yn [3],
				       wn = u->No[u->To[t][(n+0)%3 + 1]][3];

				if (wn == 0)
					wn = 1;
				else if (wm > 0)
					wn = 0;

				xn[0] = u->No[u->To[t][(n+0)%3 + 1]][1];
				xn[1] = u->No[u->To[t][(n+1)%3 + 1]][1];
				xn[2] = u->No[u->To[t][(n+2)%3 + 1]][1];
                                                         
				yn[0] = u->No[u->To[t][(n+0)%3 + 1]][2];
				yn[1] = u->No[u->To[t][(n+1)%3 + 1]][2];
				yn[2] = u->No[u->To[t][(n+2)%3 + 1]][2];

				u->A[t][m][n]  = /*wm*wn* */(ym[1] - ym[2])*(yn[1] - yn[2]);
				u->A[t][m][n] += /*wm*wn* */(xm[2] - xm[1])*(xn[2] - xn[1]);
				u->A[t][m][n] /=  (4 * T);
			}

			u->skal[t][m]  = (xm[1] - xm[0])*(ym[2] - ym[0]);
			u->skal[t][m] -= (xm[2] - xm[0])*(ym[1] - ym[0]);
			u->skal[t][m] /= (-6);
		}
	}
}

// we impose Dirichlet boundary conditions
void dirichlet (tr * u)
{
	int i, j;
	for (i = 0; i <= u->N-1; i++)
	{
		if (u->No[i][3] == 2)
		{
			// first we clean both the row and the column
			for (j = 0; j <= u->N-1; j++)
			{
				u->S[i][j] = 0;
				u->S[j][i] = 0;
			}

			// now we fix g and S[i][i]
			u->S[i][i] = 1;	// now it's not singular anymore
			u->g[i] = 0;
		}
	}
}

void calc_gS (tr * u)
{
	// at this point we have calculated A and scan
	// we only need to fill this in the matrix S and vector g
	u->S = (double **) malloc (u->N * sizeof (double *));
	u->g = (double *) calloc (u->N, sizeof (double));

	int t, m, n;
	for (t = 0; t <= u->N; t++)
		u->S[t] = (double *) calloc (u->N, sizeof (double));

	// now all is initialized so we just have to fill it
	for (t = 0; t <= u->T-1; t++)
	{
		for (m = 0; m <= 2; m++)
		{
			for (n = 0; n <= 2; n++)
				u->S[u->To[t][m+1]][u->To[t][n+1]] += u->A[t][m][n];

			u->g[u->To[t][m+1]] += u->skal[t][m];
		}
	}

	dirichlet (u);
}

void zax_fprintf (char * dat, tr * u)
{
	FILE * fout = fopen (dat, "w");
	int i;

	for (i = 0; i <= u->N-1; i++)
		fprintf (fout, "% .12e\t % .12e\t % .12e\n",
				u->v[i].x, u->v[i].y, u->v[i].val);

	fclose (fout);
}

// solve for c
void c_solver (tr * u)
{
	gsl_matrix * S = gsl_matrix_alloc (u->N, u->N);
	gsl_vector * g = gsl_vector_alloc (u->N);
	gsl_vector * c = gsl_vector_alloc (u->N);
	gsl_permutation * p = gsl_permutation_alloc (u->N);
	int signum;

	int i, j;
	for (i = 0; i <= u->N-1; i++)
	{
		gsl_vector_set (g, i, u->g[i]);
		for (j = 0; j <= u->N-1; j++)
			gsl_matrix_set (S, i, j, u->S[i][j]);
	}

	gsl_linalg_LU_decomp (S, p, &signum);
	gsl_linalg_LU_solve (S, p, g, c);
	u->c = (double *) malloc (u->N * sizeof(double));

	gsl_vector_fprintf (stdout, c, "% e");

	for (i = 0; i <= u->N-1; i++)
		u->c[i] = gsl_vector_get (c, i);

	gsl_permutation_free (p);
	gsl_vector_free (c);
	gsl_vector_free (g);
	gsl_matrix_free (S);

	u->v = (zax *) malloc (u->N * sizeof (zax));
	for (i = 0; i <= u->N-1; i++)
	{
		u->v[i].x	= u->No[i][1];
		u->v[i].y	= u->No[i][2];
		u->v[i].val	= u->c[i];
		u->v[i].attribute = u->No[i][3];
	}

	zax_fprintf ("solution.dat", u);
}

#endif

