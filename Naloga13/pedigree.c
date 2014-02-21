#include <stdio.h>
#include <stdlib.h>

int main (int argc, char ** argv)
{
	if (argc != 2)
	{
		printf ("Error! Specify the input file!\n");
		exit (EXIT_FAILURE);
	}

	size_t X, Y, i, j, k;
	FILE * fout = fopen(argv[1], "r");
	if (!fout)
	{
		printf ("Error! File does not exist.\n");
		exit (EXIT_FAILURE);
	}

	// we read some simple parameters
	fscanf (fout, "# ImageMagick pixel enumeration: %lu,%lu,%*d,gray\n",
			&X, &Y);

	// we create our array
	int ** A = (int **) malloc (X * sizeof(int *));
	for (i = 0; i <= X-1; i++)
		A[i] = (int *) malloc (Y * sizeof(int));

	int a;

	for (k = 0; k <= X*Y-1; k++)
	{
		fscanf (fout, "%lu,%lu: (%d,%*d,%*d) %*s gray(%*d,%*d,%*d)\n",
				&i, &j, &a);
		A[i][j] = a;
	}

	// now we print that bad boy out
	for (j = 0; j <= Y-1; j++)
	{
		for (i = 0; i <= X-1; i++)
			printf ("%d ", A[i][j]);

		printf ("\n");
	}
	exit (EXIT_SUCCESS);
}

