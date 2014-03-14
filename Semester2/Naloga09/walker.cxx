#include "walker.h"

Walker::Walker (int number_of_vertices)
{
	// first we allocate them
	N = number_of_vertices;
	c (N-1);
	u (N-1);
	A (N-1, N-1);
	
	// and now we set the values
	fillu ();
	fillA ();
}

Matrix2d
Walker::rotMat (int j)
{
	Vector2d t = point[j+1] - point[j];
	t.normalize ();

	Matrix2d R;
	R << t[0], -t[1],
	     t[1], t[0];

	return R;
}

double
Walker::potential (int j, Vector2d r)
{
	Matrix2d R = rotMat (j);

	Vector2d v1	= R * point[j],
		 v2	= R * point[j+1],
		 rr	= R * r;

	double x1	= rr[0] - v1[0],
	       x2	= rr[0] - v2[0],
	       y	= rr[1];

	double U = (2*x2 - 2*x1 + 2*y * (atan(x1/y) - atan(x2/y)) +
		x1*log(x1*x1 + y*y) - x2*log(x2*x2 - y*y)) / (4 * M_PI);

	return U;
}

Vector2d
Walker::velocity (int j, Vector2d r)
{
	Matrix2d R = rotMat (j);

	Vector2d v1	= R * point[j],
		 v2	= R * point[j+1],
		 rr	= R * r,
		 v;

	double x1	= rr[0] - v1[0],
	       x2	= rr[0] - v2[0],
	       y	= rr[1];

	v[0] = log((x1*x1 + y*y)/(x2*x2 + y*y)) / (4 * M_PI);
	v[1] = (atan(x1/y) - atan(x2/y)) / (2 * M_PI);

	v = R.transpose()*v;
	return v;
}

void
Walker::fillA (void)
{
	for (int i = 0; i <= N-2; i++)
	{
		for (int j = 0; j <= N-2; j++)
			A(i,j) = potential (i, (point[j] + point[j+1])/2);
	}
}

