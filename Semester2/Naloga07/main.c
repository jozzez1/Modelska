#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include "hed.h"

int main (int argc, char ** argv)
{
	int seg     = 100,
	    sgg     = 0,
	    arg;

	struct option longopts[] =
	{
		{ "seg",        required_argument,     NULL,    's' },
		{ "sgg",	no_argument,           NULL,    'p' },
		{ "help",       no_argument,           NULL,    'h' },
		{ NULL,         0,                     NULL,      0 }
	};

	while ((arg = getopt_long (argc, argv, "s:hp", longopts, NULL)) != -1)
	{
		switch (arg)
		{
			case 's':
				seg = atoi (optarg);
				break;
			case 'p':
				sgg = 1;
				break;
			case 'h':
				printf ("List of commands:\n");
				printf ("--seg=<Int>, -s<Int>\n");
				printf ("\t No. of initialsegments for triangulation.\n");
				printf ("--sgg\n");
				printf ("\t Create files for processing with the \"Triangle\" program\n");
				printf ("--help\n\t Printf this list\n");
				exit (EXIT_SUCCESS);
			default:
				printf ("Unknown command!\nTry %s --help, for list of commands\n", argv[0]);
				exit (EXIT_FAILURE);
		}
	}

	if (sgg)
	{
		PolyFile (seg);
		exit (EXIT_SUCCESS);
	}

	printf ("Initializing triangulator\n");
	tr * u = (tr *) malloc (sizeof (tr));
	printf ("Getting triangles\n");
	getTriangles (seg, u);
	printf ("Calculating A and scal\n");
	calculateAscal (u);
	printf ("Writing it into gS\n");
	calc_gS (u);
	printf ("Computing LU and c\n");
	c_solver (u);

	destroy_tr (u);

	exit (EXIT_SUCCESS);
}

