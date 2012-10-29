/* standard C libraries */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* gsl libraries */
#include <gsl/gsl_matrix.h> // matrices ...
#include <gsl/gsl_rng.h>    // random number generation
#include <gsl/gsl_blas.h>   // blas support for matrix multiplication

/* this function is always needed */
void change (int a, int b)
{
	int c = a;

	a = b;
	b = c;
}

/* my struct definition */
typedef struct
{
	int N, 		// matrix dimension -- number of nodes
	    D,		// number of connections between nodes, even number
	    S,		// random generator seed
	    diam,	// diameter of the graph
	    Imax;	// maximum allowed iteration

	gsl_rng * r;	// random number generator
	gsl_matrix * A; // current connection matrix

	double ** m,	// current connections
	       clust;	// clustering


	FILE * fout,	// output file of the clustering etc.
	     * fmatrix;	// the end-matrix
	
	char * dat,	// output filename
	     * mat;	// matrix filename
} dub;

/* allocator */
dub * alloc (void)
{
	dub * u = (dub *) malloc (sizeof (dub *));

	return u;
}

/* for debugging purposes */
void print_matrix (gsl_matrix * A)
{
	printf ("\n");
	for (int i = 0; i <= A->size1-1; i++)
	{
		for (int j = 0; j <= A->size1-1; j++)
			printf ("%d  ", (int) gsl_matrix_get(A, i, j));
		printf ("\n");
	}
}

/* to output the matrix */
void fprint_matrix (dub * u)
{
	fprintf (u->fmatrix, "\n");
	for (int i = 0; i <= u->N-1; i++)
	{
		for (int j = 0; j <= u->N-1; j++)
			fprintf (u->fmatrix, "%d  ", (int) gsl_matrix_get(u->A, u->N-1 - i, j));
		fprintf (u->fmatrix, "\n");
	}
}

/* code to set the matrix elements */
void matrix_set_init (dub * u)
{
	if (u->D%2 != 0 && u->D < 2 && u->D < u->N)
	{
		printf ("D must be a positive even number!\n");
		exit (EXIT_FAILURE);
	}
	/* we assume we start with a clean matrix */
	int i, j;
	for (i = 0; i <= u->N - 1; i++)
	{
		/* we set up the trivial connection */
		gsl_matrix_set (u->A, i, i, 1.0);

		/* now we set up all the other */
		for (j = 1; j <= u->D/2; j++)
		{
			/* we fill the matrix */
			gsl_matrix_set (u->A, (i-j+u->N)%u->N, i, 1.0);
			gsl_matrix_set (u->A, (i+j)%u->N, i, 1.0);

			/* array of all connections */
			u->m[i + (j-1)*u->N][0] = i;		// we just need the lower diagonal
			u->m[i + (j-1)*u->N][1] = (i+j)%u->N;
		}
	}
}

/* calculation of graph diameter */
void diameter_cal (dub * u)
{
	int k = 1, 	// maximum length in the graph
	    i, j,	// loop indices
	    cont = 0;	// continue flag for the 'do-while' loop

	gsl_matrix * M = gsl_matrix_alloc (u->N, u->N),
		   * R = gsl_matrix_alloc (u->N, u->N);

	/* we copy the matrix */
	gsl_matrix_memcpy (M, u->A);

	do
	{
		/* we scan for the zeroes, we only check the upper triangle, without diagonal */
		for (i = 0; i <= u->N - 1; i++)
		{
			cont = 0;
			for (j = i+1; j <= u->N - 1; j++)
			{
				if (!gsl_matrix_get(M, i, j))
				{
					cont = 1;
					break;
				}
			}
			if (cont) break;
		}

		if (cont)
		{
			gsl_blas_dsymm (CblasRight, CblasLower, 1.0, M, u->A, 0.0, R);
			gsl_matrix_memcpy (M, R);
			k++;
		}

		else break;

	} while (cont && k <= u->N-1);

	gsl_matrix_free (M);
	gsl_matrix_free (R);

	u->diam = k;
}

/* clustering coefficient calculation */
void cluster_cal (dub * u)
{
	/* 
	 * definition of the clustering coefficient
	 * C = no. of closed triplets / no. of connected triplets
	 * based on checking if neighbour of my neighbour, is my neigbour too
	 */

	int i, j, k,
	    ctriag = 0,
	    otriag = 0;

	/* we just count the lower triangle, otherwise
	 * we would count each closed triangle twice */
	for (i = 0; i <= u->N-1; i++)
	{
		for (j = i+1; j <= u->N-1; j++)
		{

			if (gsl_matrix_get (u->A, j, i))
			{
				for (k = j+1; k <= u->N-1; k++)
				{
					double p1 = gsl_matrix_get (u->A, k, i),
					       p2 = gsl_matrix_get (u->A, k, j);

					if (p1 == 1 && p2 == 1)
					{
						ctriag++;
						otriag++;
					}

					else if (p1 == 1 || p2 == 1)
						otriag++;
				}
			}
		}
	}
	u->clust = 1.0*ctriag/otriag;
}

/* initializor */
void init (dub * u, int N, int D, int S, int Imax)
{
	int i;

	u->N = N;
	u->D = D;
	u->S = S;
	u->Imax = Imax;

	if (u->S == 0)
		u->S = 5489;

	u->r = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (u->r, u->S);

	/* we prepare the current connetions' array */
	u->m = (double **) malloc (u->N*(u->D/2) * sizeof(double *));
	for (i = 0; i <= u->N*u->D/2 - 1; i++)
		u->m[i] = (double *) malloc (2 * sizeof(double));

	/* we create a matrix and initialize its elements to zero */
	u->A = gsl_matrix_calloc (u->N, u->N);

	/* we set the matrix and the connections' array*/
	matrix_set_init (u);

	/* distance of the graph */
	diameter_cal (u);

	/* clustering */
	cluster_cal (u);

	/* we open the output file for writing */
	u->dat = (char *) malloc (20);
	u->mat = (char *) malloc (20);

	sprintf (u->dat, "connectivity-%d-%d.txt", N, D);
	sprintf (u->mat, "matrix-%d-%d.txt", N, D);

	u->fout = fopen (u->dat, "w");
	u->fmatrix = fopen (u->mat, "w");
}

/* destructor */
void destroy (dub * u)
{
	gsl_rng_free (u->r);
	gsl_matrix_free (u->A);
	
	fclose (u->fout);
	fclose (u->fmatrix);

	int i;
	for (i = 0; i <= u->N*u->D/2 - 1; i++)
		free (u->m[i]);

	free (u->m);

	free (u->dat);
	free (u->mat);

	free (u);
}

/* propagator */

/*
 * the point of the propagator is this: it scans over the
 * existing connections, delets one, and creates it somewhere
 * else and it does so in a different but simple ways
 * 
 */

void step (dub * u)
{
	/* first we pull a random number from our "hat" */
	int rand, // random connection index
	    i, j; // nodes of the connection

	do
	{
		rand = gsl_rng_uniform_int (u->r, u->N*u->D/2);
		i = u->m [rand][0]; // connection goes from i to j
		j = u->m [rand][1]; // ...
	} while (i == j);

	/* we erase the value in the matrix :D ... trololololo... */
	if (gsl_matrix_get (u->A, i, j) != 1)
	{
		printf ("Error!\nA non-unity element!\n");
		printf ("Defect matrix construction.\n");
		exit (EXIT_FAILURE);
	}

	else
	{
		/* we erase that 1 */
		gsl_matrix_set (u->A, i, j, 0);
		gsl_matrix_set (u->A, j, i, 0);

		/* if the probability for finding a hole randomly is larger than 0.5
		 * then let's do it randomly */
		if (u->N >= 2*u->D)
		{
			do
			{
				i = gsl_rng_uniform_int (u->r, u->N);
				j = gsl_rng_uniform_int (u->r, u->N);
			} while (gsl_matrix_get(u->A, i, j) || i == j);
		}

		/* if doing randomly doesn't pay off, we will do it like men do!
		 * with brute force! */
		else
		{
			do
			{
				int counter = 1;
				do
				{
					i = (i+1)%u->N;
					j = (j+1)%u->N;
					counter++;
				} while (counter != u->N && gsl_matrix_get(u->A, i, j));

				if (!gsl_matrix_get(u->A, i, j) && i != j)
					break;
				else
					j++;
				if (i == j)
					j++;
			} while (gsl_matrix_get(u->A, i, j));
		}

		/* put that 1 in the matrix */
		gsl_matrix_set (u->A, i, j, 1);
		gsl_matrix_set (u->A, j, i, 1);

		u->m[rand][0] = i;
		u->m[rand][1] = j;

		/* we fix the parameters */
		diameter_cal (u);
		cluster_cal (u);
	}
}

void dump (dub * u, int i)
{
	fprintf (u->fout, "% d\t% d\t% .3lf\n", i, u->diam, u->clust);
}

double norm (int diam, double clust)
{
	double re = sqrt(diam*diam + clust*clust);

	return re;
}

void propagate (dub * u)
{
	dump (u, 0);

	int diam_old = 0,
	    diam_new = u->diam;

	double clust_old = 0,
	       clust_new = u->clust;

	int i = 0;
	while (u->clust >= 0.01 && i <= u->Imax)
	{
		diam_old = diam_new;
		clust_old = clust_new;

		step (u);
		i++;
		dump (u,i);
		if (i%500 == 0)
		{
			printf("i = %d\n", i);
			clust_new = u->clust;
			diam_new = u->diam;

			if (clust_new > clust_old)
				break;

			else if (norm (diam_old - diam_new, clust_old - clust_new) <= 0.001)
				break;
		}
	}
	fprint_matrix (u);
}

int main (int argc, char ** argv)
{
	if (argc > 4 || argc < 3)
	{
		printf ("Wrong program usage!\n");
		printf ("Correct usage:\n");
		printf ("%s <N> <D> [<S>]\n", argv[0]);

		exit (EXIT_FAILURE);
	}

	else
	{
		int N = atoi (argv[1]),
		    D = atoi (argv[2]),
		    I = 10000,
		    S = 0;
		if (argc == 4)
			S = atoi (argv[3]);

		dub * u = alloc ();
		init (u, N, D, S, I);
		propagate (u);


		printf ("Finished!\n");
		printf ("%s created!\n%s created!\n", u->dat, u->mat);

		destroy (u);
		exit (EXIT_SUCCESS);
	}
}
