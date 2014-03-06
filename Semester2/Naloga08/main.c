#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include "hed.h"

int main (int argc, char ** argv)
{
	double area	= 0.01,
	       angle	= 20;

	int seg		= 100,
	    w_eig	= 0,
	    arg;

	struct option longopts[] =
	{
		{ "seg",        required_argument,     NULL,    's' },
		{ "area",       required_argument,     NULL,    'a' },
		{ "angle",      required_argument,     NULL,    'q' },
		{ "solve",      no_argument,           NULL,    'y' },
		{ "help",       no_argument,           NULL,    'h' },
		{ NULL,         0,                     NULL,      0 }
	};

	while ((arg = getopt_long (argc, argv, "s:a:q:hy", longopts, NULL)) != -1)
	{
		switch (arg)
		{
			case 's':
				seg = atoi (optarg);
				break;
			case 'a':
				area = atof (optarg);
				break;
			case 'q':
				angle = atof (optarg);
				break;
			case 'y':
				w_eig = 1;
				break;
			case 'h':
				printf ("List of commands:\n");
				printf ("--seg=<Int>, -s<Int>\n");
				printf ("\t No. of initialsegments for triangulation.\n\n");
				printf ("--area=<Double>, -a<Double>\n");
				printf ("\t Maximum triangle area.\n\n");
				printf ("--angle=<Double>, -q<Double>\n");
				printf ("\t Minimum triangle qngle.\n\n");
				printf ("--help\n\t Printf this list\n");
				exit (EXIT_SUCCESS);
			default:
				printf ("Unknown command!\nTry %s --help, for list of commands\n", argv[0]);
				exit (EXIT_FAILURE);
		}
	}

	tr * u = (tr *) malloc (sizeof (tr));
	u->w_eig = w_eig;
	solve (seg, u, angle, area);

	destroy_tr (u);

	exit (EXIT_SUCCESS);
}

