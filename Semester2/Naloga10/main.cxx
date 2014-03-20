#include <iostream>
#include <fstream>
#include <string>
#include <mgl2/data.h>
#include <mgl2/mgl.h>

int main (int argc, char ** argv)
{
	if (argc != 3)
	{
		std::cout << "Wrong program usage! Try:" 	<< std::endl;
		std::cout << argv[0] << " <filename> <N>" 	<< std::endl;
		exit (EXIT_FAILURE);
	}

	std::string filename (argv[1]);
	int N = atoi (argv[2]);

	double * x = new double[N*N],
		   * y = new double[N*N],
		   * z = new double[N*N];

	std::ifstream fin (filename.c_str(), std::ifstream::in);
	if (!fin)
	{
		std::cerr << "Error. File not found." << std::endl;
		exit (EXIT_FAILURE);
	}

	for (int i = 0; !fin.eof(); i++)
	{
		fin >> x[i] >> y[i] >> z[i];
	}
	fin.close ();

	std::cout << x[2] << y[2] << z[2] <<  std::endl;

	exit (EXIT_SUCCESS);
}
