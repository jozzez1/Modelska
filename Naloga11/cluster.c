#include <stdio.h>
#include <stdlib.h>

int main (int argc, char ** argv)
{
	if (argc != 2)
		exit (EXIT_FAILURE);

	// OK let's open the file
	int N = atoi (argv[1]),
	    i = 0,
	    j = 0,
	    k = 0;

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

	// so everything seems to be in order, let's calculate
	// the clustering with brute force
	double C = 0;
	for (i = 0; i <= N-1; i++)
	{
		for (j = 0; j <= N-1; j++)
		{
			if (j == i) continue;
			for (k = 0; k <= N-1; k++)
			{
				if ((k == i) || (k == j)) continue;

				if ((M[j][i] == 1) && (M[k][i] == 1) && (M[k][j] == 1))
					C++;
			}
		}
	}

	// now we have to free the allocated space
	for (i = 0; i <= N-1; i++)
		free (M[i]);
	free (M);

	C /= (N * (N-1) * (N-2));
	printf ("%e\n", C);

	exit (EXIT_SUCCESS);
}

