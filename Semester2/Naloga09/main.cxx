#include <iostream>
#include <getopt.h>
#include "walker.h"

int main (int argc, char ** argv)
{
	int N     = 200,
	    mode  = 0,
	    nx    = 200,
	    ny    = 200,
	    arg;
	double a  = 0,
	       b  = 0,
	       xm = -5,
	       xM = +5,
	       ym = -5,
	       yM = +5,
	       ui = 1;

	struct option longopts[] =
	{
		{ "Number", required_argument,       NULL,         'N' },
		{ "mode",   required_argument,       NULL,         'm' },
		{ "p1",     required_argument,       NULL,         'a' },
		{ "p2",     required_argument,       NULL,         'b' },
		{ "nx",     required_argument,       NULL,         'I' },
		{ "ny",     required_argument,       NULL,         'J' },
		{ "xmin",   required_argument,       NULL,         'x' },
		{ "xmax",   required_argument,       NULL,         'X' },
		{ "ymin",   required_argument,       NULL,         'y' },
		{ "ymax",   required_argument,       NULL,         'Y' },
		{ "uinf",   required_argument,       NULL,         'u' },
		{ "help",   no_argument,             NULL,         'h' },
		{0, 0, 0, 0}
	};

	while ((arg = getopt_long (argc, argv, "N:m:a:b:I:J:x:X:y:Y:u:h", longopts, NULL)) != -1)
	{
		switch (arg)
		{
			case 'N': N	= atoi (optarg); break;
			case 'm': mode	= atoi (optarg); break;
			case 'a': a	= atof (optarg); break;
			case 'b': b	= atof (optarg); break;
			case 'I': nx	= atoi (optarg); break;
			case 'J': ny	= atoi (optarg); break;
			case 'x': xm	= atof (optarg); break;
			case 'X': xM	= atof (optarg); break;
			case 'y': ym	= atof (optarg); break;
			case 'Y': yM	= atof (optarg); break;
			case 'u': ui	= atof (optarg); break;
			case 'h':
				  std::cout << "-N, --Number:  set number of vertices" 		<< std::endl;
				  std::cout << "-m, --mode:    set mode [0-3]" 			<< std::endl;
				  std::cout << "-a, --p1:      set the 1st mode parameter" 	<< std::endl;
				  std::cout << "-b, --p2:      set the 2nd mode parameter" 	<< std::endl;
				  std::cout << "-I, --nx:      set the number of 'x' points"	<< std::endl;
				  std::cout << "-J, --ny:      set the number of 'y' points"	<< std::endl;
				  std::cout << "-x, --xmin:    set the minimum 'x'"		<< std::endl;
				  std::cout << "-X, --xmax:    set the maximum 'x'"		<< std::endl;
				  std::cout << "-y, --ymin:    set the minimum 'y'"		<< std::endl;
				  std::cout << "-Y, --ymax:    set the maximum 'y'"		<< std::endl;
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

	Walker * u = new Walker (N, mode, a, b,
			nx, ny, xm, xM, ym, yM,
			ui);
	u->solve (mode);

	u->~Walker();

	exit (EXIT_SUCCESS);
}

