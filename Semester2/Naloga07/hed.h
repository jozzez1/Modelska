// hed.h
///////////
#ifndef __HEADER_MOD207
#define __HEADER_MOD207

#include <math.h>
#include <stdlib.h>

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
	long int ** To;
	
	double *** A,
	       ** skal,
	       ** S,
	       * g;
} tr;

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

	if (segments)
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
	fscanf (fnode, "%d %*d %*d\n", &N);

	u->N = N;
	u->No = (double **) malloc (N * sizeof (double *));
	for (i = 0; i <= N-1; i++)
	{
		u->No[i] = (double *) malloc (4 * sizeof (double));
		fscanf (fnode, "%lf %lf %lf %lf\n",
				&u->No[i][0],	// vertice index
				&u->No[i][1],	// "x" component
				&u->No[i][2],	// "y" component
				&u->No[i][3]);	// if we are on the edge or not
	}
	fclose (fnode);

	// an now the triangles
	fscanf (fele, "%d %*d %*d\n", &N);

	u->T = N;
	u->To = (long int **) malloc (N * sizeof (long int *));
	for (i = 0; i <= N-1; i++)
	{
		u->To[i] = (long int *) malloc (4 * sizeof (long int));
		fscanf (fele, "%lu %lu %lu %lu\n",
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
void calculateA (tr * u)
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
			       ym [3];

			xm[0] = u->No[u->To[t][(m+0)%3]][1];
			xm[1] = u->No[u->To[t][(m+1)%3]][1];
			xm[2] = u->No[u->To[t][(m+2)%3]][1];

			ym[0] = u->No[u->To[t][(m+0)%3]][2];
			ym[1] = u->No[u->To[t][(m+1)%3]][2];
			ym[2] = u->No[u->To[t][(m+2)%3]][2];

			for (n = 0; n <= 2; n++)
			{
				double xn [3],
				       yn [3];

				xn[0] = u->No[u->To[t][(n+0)%3]][1];
				xn[1] = u->No[u->To[t][(n+1)%3]][1];
				xn[2] = u->No[u->To[t][(n+2)%3]][1];
                                                         
				yn[0] = u->No[u->To[t][(n+0)%3]][2];
				yn[1] = u->No[u->To[t][(n+1)%3]][2];
				yn[2] = u->No[u->To[t][(n+2)%3]][2];

				u->A[t][m][n]  = (ym[1] - ym[2])*(yn[1] - yn[2]);
				u->A[t][m][n] += (xm[2] - xm[1])*(xn[2] - xn[1]);
				u->A[t][m][n] /=  (4 * T);
			}

			u->skal[t][m]  = (xm[1] - xm[0])*(ym[2] - ym[0]);
			u->skal[t][m] -= (xm[2] - xm[0])*(ym[1] - ym[0]);
			u->skal[t][m] /= (-6);
		}
	}
}

#endif

