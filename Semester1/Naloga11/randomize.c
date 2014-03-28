#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double clustering (int ** M, const int N)
{
	int i, j, k;
	double C = 0,
	       T, K;
	for (i = 0; i <= N-1; i++)
	{
		T = 0;
		K = 0;
		// count the triangles on 'i'
		for (j = 0; j <= N-1; j++)
		{
			if (j == i) continue;
			for (k = 0; k <= N-1; k++)
			{
				if ((k == i) || (k == j)) continue;

				if ((M[j][i] == 1) && (M[k][i] == 1) && (M[k][j] == 1))
					T++;
			}

			K += M[j][i];
		}
		K = K*(K-1)/2;
		T /= 2;

		if (K != 0)
			C += (T/K)/N;
	}
	return C;
}

void product (int ** C, const int ** A, const int ** B, const int N)
{
	// C = A * B
	int i, j, k;
	for (i = 0; i <= N-1; i++)
	{
		for (j = 0; j <= N-1; j++)
		{
			C[i][j] = 0;
			for (k = 0; k <= N-1; k++)
				C[i][j] += A[i][k] * B[k][j];
		}
	}
}

void equals (int ** A, const int ** B, const int N)
{
	int i, j;
	for (i = 0; i <= N-1; i++)
	{
		for (j = 0; j <= N-1; j++)
			A[i][j] = B[i][j];
	}
}

void unit (int ** A, const int N)
{
	int i, j;
	for (i = 0; i <= N-1; i++)
	{
		for (j = 0; j <= N-1; j++)
			A[i][j] = 0;
	}

	for (i = 0; i <= N-1; i++)
		A[i][i] = 1;
}

int look4zero (int ** M, const int N)
{
	int i, j,
	    found = 0;
	for (i = 0; i <= N-1; i++)
	{
		for (j = 0; j <= N-1; j++)
		{
			if (M[i][j] == 0)
				found = 1;
		}
	}
	return found;
}

int sum (int ** A, const int j, const int N)
{
	int S = 0,
	    i;
	for (i = 0; i <= N-1; i++)
	{
		S += A[i][j];
	}
	return S;
}

void randcon (int ** M, const int N)
{
	size_t i, j, k, S,
	    kx = 0,
	    ky = 0;
	do
	{
		j = random () % N;
		S = sum (M, j, N);
	} while (S == N || S == 1);

	long int * x = (long int *) malloc ((kx + 1) * sizeof (long int)),
	    * y = (long int *) malloc ((ky + 1) * sizeof (long int));
	
	for (i = 0; i <= N-1; i++)
	{
		if (M[i][j] == 1 && i != j)
		{
			x[kx] = i;
			kx++;
			x = (long int *) realloc (x, (kx + 1) * sizeof (long int));
		}
		else if (M[i][j] == 0)
		{
			y[ky] = i;
			ky++;
			y = (long int *) realloc (y, (ky + 1) * sizeof (long int));
		}
	}

	i = x [random() % kx];
	k = y [random() % ky];

	if (x) free (x);
	if (y) free (y);

	M [i][j] = 0;
	M [j][i] = 0;

	M [k][j] = 1;
	M [j][k] = 1;
}

int radiusz (int ** M, const int N)
{
	int r = 0,
	    i;

	int ** A = (int **) malloc (N * sizeof (int *)),
	    ** B = (int **) malloc (N * sizeof (int *));
	for (i = 0; i <= N-1; i++)
	{
		A[i] = (int *) malloc (N * sizeof (int));
		B[i] = (int *) malloc (N * sizeof (int));
	}

	unit (A, N);
	
	do
	{
		equals (B, (const int **) A, N);
		product (A, (const int **) B, (const int **) M, N);
		r++;
	} while (look4zero (A, N));

	for (i = 0; i <= N-1; i++)
	{
		free (A[i]);
		free (B[i]);
	}
	free (A);
	free (B);

	return r;
}

int main (int argc, char ** argv)
{
	if (argc != 4)
		exit (EXIT_FAILURE);

	// OK let's open the file
	const int N = atoi (argv[1]),
	          p = atoi (argv[2]);

	double percent = atof (argv[3]);

	int i = 0,
	    j = 0;

	char c, dat [40];
	sprintf (dat, "matrix-%d.txt", N);

	FILE * fin = fopen (dat, "r");
	if (!fin)
		exit (EXIT_FAILURE);

	// everything seems to exist .. let's create the matrix
	int ** M = (int **) malloc (N * sizeof (int *));
	for (i = 0; i <= N-1; i++)
		M[i] = (int *) malloc (N * sizeof (int));

	// and now read the elements from the file
	i = 0, j = 0;
	do
	{
		c = fgetc (fin);
		if (c == '1' || c == '0')
		{
			M [i][j] = atoi (&c);
			j++;

		}
		if (c == '\n')
		{
			i++;
			j = 0;
		}
	} while (c != EOF);
	fclose (fin);

	// so everything seems to be in order, let's randomize and get the plots
	srandom (time(NULL));

	
	int D = 0;
	for (i = 0; i <= N-1; i++)
		D += sum (M, i, N);
	D = D/N-1;

	FILE * fout;
	char dat1 [40];
	sprintf (dat1, "progress-N%d-D%d.txt",N, D);
	if (p) fout = fopen (dat1, "w");

	int Imax = N*D*percent/2;
	printf ("Imax = %d\n", Imax);
	for (i = 0; i <= Imax-1; i++)
	{
		randcon (M, N);
//		printf ("%d/%d\n", i, Imax-1);

		if (p)
		{
			int r = radiusz (M, N);
			double C = clustering (M, N);
		
			fprintf (fout, "%d\t %d\t %lf\n",
					i, r, C);
		}
	}

	if (p) fclose (fout);

	// the final output
	fin = fopen (dat, "w+");

	for (i = 0; i <= N-1; i++)
	{
		for (j = 0; j <= N-1; j++)
			fprintf (fin, "%d ", M[i][j]);
		fprintf (fin, "\n");
	}

	fclose (fin);

	// now we have to free the allocated space
	for (i = 0; i <= N-1; i++)
		free (M[i]);
	free (M);

	exit (EXIT_SUCCESS);
}

