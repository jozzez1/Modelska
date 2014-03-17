#include "profiles.h"

// just a simple line: $x \in [0,1],\ y = 0$
std::vector<Eigen::Vector2d>
simple (int N, void * params)
{
	// we don't really need params :P
	double y = *((double *) params);
	y = 0;
	std::vector<Eigen::Vector2d> points;
	for (int i = 0; i <= N-1; i++)
		points.push_back (Eigen::Vector2d (i*1.0/(N-1), y));

	return points;
}

// elipsoidal profile
std::vector<Eigen::Vector2d>
elipsoid (int N, void * params)
{
	double b = *((double *) params);
	std::vector<Eigen::Vector2d> points;
	for (int i = 0; i <= N-1; i++)
	{
		double x = cos (2*M_PI*i/(N-1)),
		       y = b * sin (2*M_PI*i/(N-1));

		points.push_back (Eigen::Vector2d (x, y));
	}

	return points;
}

// fish profile
std::vector<Eigen::Vector2d>
fishpr (int N, void * params)
{
	assert (N & 1); // test if N is even

	double t = *((double *) params);
	std::vector<Eigen::Vector2d> points;
	for (int i = 0; i <= (N+1)/2 - 1; i++)
	{
		double x = 1 - (1.0*i)/((N+1)/2 - 1),
		       y = 1.457122 * sqrt(x) 
			       - 0.624424*x
			       - 1.727016*x*x
			       + 1.384087*pow(x, 3)
			       - 0.489769*pow(x, 4);

		y *= t/50;
		points.push_back (Eigen::Vector2d (x, y));
	}

	for (int i = 1; i <= (N+1)/2 - 1; i++)
	{
		double x = (1.0*i)/((N+1)/2 - 1),
		       y = -1.457122 * sqrt(x) 
			       + 0.624424*x
			       + 1.727016*x*x
			       - 1.384087*pow(x, 3)
			       + 0.489769*pow(x, 4);

		y *= t/50;
		points.push_back (Eigen::Vector2d (x, y));
	}

	points.push_back (points[0]);

	return points;
}

// Zukovski profile
std::vector<Eigen::Vector2d>
zukovski (int N, void * params)
{
	// we transform a unit circle into something else 
	double A = ((double *) params) [0],
	       B = ((double *) params) [1];

	std::vector<Eigen::Vector2d> points;
	for (int i = 0; i <= N-1; i++)
	{
		double x  = A + 0.9*cos (2 * M_PI * i/(N-1)),
		       y  = B + 0.9*sin (2 * M_PI * i/(N-1)),
		       n  = x*x + y*y,
		       rx = 0.5*(x + x/n),
		       ry = 0.5*(y - y/n);

		points.push_back (Eigen::Vector2d (rx, ry));
	}

	return points;
}
