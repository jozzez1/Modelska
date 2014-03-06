// hed.h
///////////
#ifndef __HEADER_MOD207
#define __HEADER_MOD207

#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_eigen.h>

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
	       N,	// number of vertices
	       n;	// number of vertices that are not on the boundary -- Dirichlet

	zax ** p;	// height at each vertex

	double ** No;	// vertex matrix
	int ** To,	// triangle matrix
	    * no;	// vector of indices of non-edge vertices -- Dirichlet
	
	double *** A,	// matrix A -- look Sirca
	       ** S,	// matrix S
	       ** B,	// matrix B
	       * e;	// eigenvalues
	
	double a, q;

	int w_eig;
} tr;

void destroy_tr (tr * u)
{
	return;
	int i, j;
	// free the vector
	if (u->e) free (u->e);
	if (u->no) free (u->no);

	// now the matrices
	for (i = 0; i <= u->N-1; i++)
	{
		if (u->No[i]) free (u->No[i]);
		if (u->S[i]) free (u->S[i]);
		if (u->B[i]) free (u->B[i]);
	}
	if (u->No) free (u->No);
	if (u->S) free (u->S);
	if (u->B) free (u->B);

	for (i = 0; i <= 11; i++)
		if (u->p[i]) free (u->p[i]);
	if (u->p) free (u->p);

	for (i = 0; i <= u->T-1; i++)
	{
		if (u->To[i]) free (u->To[i]);
		for (j = 0; j <= 2; j++)
			if (u->A[i][j]) free (u->A[i][j]);
		if (u->A[i]) free (u->A[i]);
	}
	if (u->To) free (u->To);
	if (u->A) free (u->A);

	// and lastly free the pointer
	if (u) free (u);
}

// creates the "*.poly" file
void PolyFile (int segments)
{
	FILE * fout = fopen ("file.poly", "w");
	int N, i;
	if (segments)
		N = segments;
	else
		N = 6;

	// create edges
	zax * edge = (zax *) malloc (N * sizeof (zax));

	if (!segments)
	{
		edge [0].x =  0.0; edge[0].y =  0.0; edge[0].val = 0; edge[0].attribute = 2;
		edge [1].x = 15.0; edge[1].y =  0.0; edge[1].val = 1; edge[1].attribute = 2;
		edge [2].x = 15.0; edge[2].y = 15.0; edge[2].val = 2; edge[2].attribute = 2;
		edge [3].x = 10.0; edge[3].y = 10.0; edge[3].val = 3; edge[3].attribute = 2;
		edge [4].x =  5.0; edge[4].y = 10.0; edge[4].val = 4; edge[4].attribute = 2;
		edge [5].x =  0.0; edge[5].y = 15.0; edge[5].val = 5; edge[5].attribute = 2;
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
		fprintf (fout, "%d\t %.16e\t %.16e\t %d\n",
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
	FILE * fnode = fopen ("file.1.node", "r"),
	     * fele  = fopen ("file.1.ele", "r");

	if (!fnode || !fele)
	{
		printf ("No input. Exiting.\n");
		exit (EXIT_FAILURE);
	}

	// let's do the nodes first
	fscanf (fnode, "%d %*d %*d %*d\n", &N);

	u->N = N;
	u->n = 0;
	u->No = (double **) malloc (N * sizeof (double *));
	u->no = (int *) malloc (sizeof (int));
	for (i = 0; i <= N-1; i++)
	{
		u->No[i] = (double *) malloc (4 * sizeof (double));
		fscanf (fnode, "%lf\t%lf\t%lf\t%lf\n",
				&u->No[i][0],	// vertice index
				&u->No[i][1],	// "x" component
				&u->No[i][2],	// "y" component
				&u->No[i][3]);	// if we are on the edge or not

		if (u->No[i][3] == 0)
		{
			u->n++;
			u->no = (int *) realloc (u->no, u->n * sizeof (int));
			u->no[u->n-1] = u->No[i][0];
		}
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

void Triangulate (int segments, tr * u, double q, double a)
{
	PolyFile (segments);
	// create the command stuff
	char command [60],
	     a_flag [14],
	     q_flag [14];
	
	if (a) sprintf (a_flag, "-a%.8lf", a);
	if (q) sprintf (q_flag, "-q%.8lf", q);

	sprintf (command, "triangle -Dp %s %s file.poly",
			a_flag, q_flag);

	// and here we execute it ..
	system (command);
	printf ("%s\n", command);

	getTriangles (segments, u);

	system ("showme file.1.ele");
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
void calculateA (tr * u)
{
	int t, m, n;

	u->A = (double ***) malloc (u->T * sizeof (double **));
	for (t = 0; t <= u->T-1; t++)
	{
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
		}
	}
}

void calc_SB (tr * u)
{
	// at this point we have calculated A and scan
	// we only need to fill this in the matrix S and vector g
	u->S = (double **) malloc (u->N * sizeof (double *));
	u->B = (double **) malloc (u->N * sizeof (double *));

	int t, m, n;
	for (t = 0; t <= u->N; t++)
	{
		u->S[t] = (double *) calloc (u->N, sizeof (double));
		u->B[t] = (double *) calloc (u->N, sizeof (double));
	}

	// now all is initialized so we just have to fill it
	for (t = 0; t <= u->T-1; t++)
	{
		double T = surf (u, t);
		for (m = 0; m <= 2; m++)
		{
			for (n = 0; n <= 2; n++)
			{
				u->S[u->To[t][m+1]][u->To[t][n+1]] += u->A[t][m][n];
				if (n == m) u->B[u->To[t][m+1]][u->To[t][n+1]] += T/6;
				else u->B[u->To[t][m+1]][u->To[t][n+1]] += T/12;
			}
		}
	}
}

void zax_fprintf (char * basename, tr * u)
{
	char * dat = (char *) malloc (40 * sizeof (char));
	FILE * fout;
	int i, j;

	for (j = 0; j <= 11; j++)
	{
		sprintf (dat, "%s-%d.dat",
				basename, j);

		fout = fopen (dat, "w");
		fprintf (fout, "# E = %lf\n", u->e[j]);

		for (i = 0; i <= u->N-1; i++)
			fprintf (fout, "% .12e\t% .12e\t% .12e\n",
					u->p[j][i].x, u->p[j][i].y, u->p[j][i].val);
	}
	if (dat) free (dat);
	if (fout) fclose (fout);
}

// solve for c -- this time we have generalized eigenvalue problem
void eig_solver (tr * u)
{
	gsl_matrix * A = gsl_matrix_alloc (u->n, u->n);
	gsl_matrix * B = gsl_matrix_alloc (u->n, u->n);
	gsl_matrix * v = gsl_matrix_alloc (u->n, u->n);
	gsl_vector * e = gsl_vector_alloc (u->n);

	int i, j;
	for (i = 0; i <= u->n-1; i++)
	{
		for (j = 0; j <= u->n-1; j++)
		{
			gsl_matrix_set (A, i, j, u->S[u->no[i]][u->no[j]]);
			gsl_matrix_set (B, i, j, u->B[u->no[i]][u->no[j]]);
		}
	}

	// here comes the eigenvalue problem ...
	gsl_eigen_gensymmv_workspace * wrk = gsl_eigen_gensymmv_alloc (u->n);
	gsl_eigen_gensymmv (A, B, e, v, wrk);
	gsl_eigen_gensymmv_free (wrk);
	gsl_matrix_free (A);
	gsl_matrix_free (B);
	gsl_eigen_gensymmv_sort (e, v, GSL_EIGEN_SORT_VAL_ASC);

	u->p = (zax **) malloc (12 * sizeof (zax *));
	u->e = (double *) malloc (12 * sizeof (double));
	for (i = 0; i <= 11; i++)
	{
		u->p[i] = (zax *) malloc (u->N * sizeof (zax));
		u->e[i] = gsl_vector_get (e, i);

		for (j = 0; j <= u->N-1; j++)
		{
			u->p[i][j].x	= u->No[j][1];
			u->p[i][j].y	= u->No[j][2];
			u->p[i][j].val  = 0;
			u->p[i][j].attribute = u->No[j][3];
		}
		for (j = 0; j <= u->n-1; j++)
			u->p[i][u->no[j]].val	= gsl_matrix_get (v, j, i);
	}

	gsl_matrix_free (v);
	gsl_vector_free (e);
}

void solve (int seg, tr * u, double q, double a)
{
	Triangulate (seg, u, q, a);
	calculateA (u);
	calc_SB (u);
	if (u->w_eig)
	{
		eig_solver (u);
		zax_fprintf ("mode", u);
	}

//	system ("rm -rf file.poly file.1.*");
}

#endif

