#include <iostream>
#include <getopt.h>
#include "walker.h"

int main (int argc, char ** argv)
{
	int N     = 200,
	    mode  = 0,
	    arg;
	double a  = 0,
	       b  = 0;

	struct option longopts[] =
	{
		{ "Number", required_argument,       NULL,         'N' },
		{ "mode",   required_argument,       NULL,         'm' },
		{ "p1",     required_argument,       NULL,         'a' },
		{ "p2",     required_argument,       NULL,         'b' },
		{ "help",   no_argument,             NULL,         'h' },
		{0, 0, 0, 0}
	};

	while ((arg = getopt_long (argc, argv, "N:m:a:b:h", longopts, NULL)) != -1)
	{
		switch (arg)
		{
			case 'N': N	= atoi (optarg); break;
			case 'm': mode	= atoi (optarg); break;
			case 'a': a	= atof (optarg); break;
			case 'b': b	= atof (optarg); break;
			case 'h':
				  std::cout << "-N, --Number:  set number of vertices" 		<< std::endl;
				  std::cout << "-m, --mode:    set mode [0-3]" 			<< std::endl;
				  std::cout << "-a, --p1:      set the 1st mode parameter" 	<< std::endl;
				  std::cout << "-b, --p2:      set the 2nd mode parameter" 	<< std::endl;
				  std::cout << "-h, --help:    print this list" 		<< std::endl;
				  exit (EXIT_SUCCESS);
			default:
				  std::cout << "Wrong usage!" 		<< std::endl;
				  std::cout << "Try" 			<< std::endl;
				  std::cout << argv[0] << " -h" 	<< std::endl;
				  std::cout << "to see options." 	<< std::endl;
				  exit (EXIT_FAILURE);
		}
	}

	Walker * u = new Walker (N, mode, a, b);
	u->solve (mode);

	u->~Walker();

	exit (EXIT_SUCCESS);
}

