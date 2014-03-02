// hed.h
///////////
#ifndef __HEADER_MOD207
#define __HEADER_MOD207

#include <math.h>

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

#endif
