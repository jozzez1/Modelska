#include <iostream>
#include "walker.h"

int main (void)
{
	int N     = 200,
	    mode  = 0;
	double a  = 0,
	       b  = 0;

	Walker * u = new Walker (N, mode, a, b);
	u->solve4c ();
	u->print_solution ();

	u->~Walker();

	return 0;
}

