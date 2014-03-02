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
		{ "help",       no_argument,           NULL,    'h' },
		{ NULL,         0,                     NULL,      0 }
	};

	while ((arg = getopt_long (argc, argv, "s:h", longopts, NULL)) != -1)
	{
		switch (arg)
		{
			case 's':
				seg = atoi (optarg);
				break;
			case 'h':
				printf ("List of commands:\n");
				printf ("--seg=<Int>, -s<Int>\n");
				printf ("\t No. of initialsegments for triangulation.\n");
				printf ("--help\n\t Printf this list\n");
				exit (EXIT_SUCCESS);
			default:
				printf ("Unknown command!\nTry %s --help, for list of commands\n", argv[0]);
				exit (EXIT_FAILURE);
		}
	}

	PolyFile (seg);
	exit (EXIT_SUCCESS);
}

