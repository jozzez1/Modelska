#include "walker.h"
#include "profiles.h"

// constructor
Walker::Walker (int number, int mode, double p1, double p2,
		int Nx, int Ny, double Xmin, double Xmax, double Ymin, double Ymax)
{
	// first we allocate them
	N = number;
	c.resize (N-1);
	u.resize (N-1);
	A.resize (N-1, N-1);

	nx	= Nx;
	ny	= Ny;
	xmin	= Xmin;
	xmax	= Xmax;
	ymin	= Ymin;
	ymax	= Ymax;

	switch (mode)
	{
		case 0: point_create = &simple;   break;
		case 1: point_create = &elipsoid; break;
		case 2: point_create = &fishpr;   break;
		case 3: point_create = &zukovski; break;
		default:
			std::cout << "Wrong mode." << std::endl;
			exit (EXIT_FAILURE);
	}

	init_points (N, mode, p1, p2);
	
	// and now we set the values
	fillu ();
	fillA ();
}

// destructor
Walker::~Walker(void)
{

}

void
Walker::init_points (int number, int mode, double p1, double p2)
{
	if (mode == 3)
	{
		double * params = new double [2];
		params[0] = p1;
		params[1] = p2;

		point = point_create (number, params);

		delete [] params;
	}
	else
		point = point_create (number, &p1);
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

	double U = (2*x2 - 2*x1 + 2*y * (atan(x1/y) - atan(x2/y))
			+ x1*log(x1*x1 + y*y) - x2*log(x2*x2 + y*y))/(4 * M_PI);

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

void
Walker::plot_chr (void)
{
	mglData x (N-1),
		y (N-1);

	y.Set (c.data(), N-1);
	for (int i = 0; i <= N-2; i++)
		x.a [i] = ((point[i+1] + point[i])/2)[0];

	// plot it with MathGL
	mglGraph gr;
	gr.Title ("Porazdelitev naboja, \\sigma(x)");
	gr.SetRange ('x', 0, 1);
	gr.SetRange ('y', 0, 120);
	gr.Label ('y', "\\sigma(x)");
	gr.Label ('x', "x");
	gr.Axis ();
	gr.Grid ("xy", "B;");
//	gr.Box ();
	gr.Plot (x, (-1)*y, "r");
	gr.WriteEPS ("potential.eps", "My test plot");
}

void
Walker::plot_pot (void)
{
	MatrixXd P (nx * ny, 3);

	double hx = (xmax - xmin)/nx,
	       hy = (ymax - ymin)/ny;

	for (int i = 0; i <= nx-1; i++)
	{
		for (int j = 0; j <= ny-1; j++)
		{
			P (j + ny*i, 0) = hx*i + xmin;
			P (j + ny*i, 1) = hy*j + ymin;
			P (j + ny*i, 2) = 0;
			for (int k = 0; k <= N-2; k++)
				P (j + ny*i, 2) += potential (k, P.block(j + ny*i, 0, 1, 2).transpose());
		}
	}

	mglGraph gr;

	mglData x (nx, ny),
		y (nx, ny),
		z (nx, ny);

	x.Set (P.block(0, 0, nx*ny, 1).data(), nx, ny);
	y.Set (P.block(0, 1, nx*ny, 1).data(), nx, ny);
	z.Set (P.block(0, 2, nx*ny, 1).data(), nx, ny);

	gr.SetSize (800, 800);
	gr.Title ("Potencial - U(x,y)");
	gr.SetDefScheme ("wyqrRk");
	gr.SetRanges (x, y);
	gr.SetRange ('c', z);
	gr.Colorbar ("wyqrRk^");
	gr.Axis ();
	gr.Label ('x', "x");
	gr.Label ('y', "y");
	gr.Label ('c', "U(x,y)");
	gr.Box ();

	gr.Dens (x, y, z);
	gr.Aspect (1,1);
	gr.WritePNG ("Potencial.png");
}

void
Walker::solve (int mode)
{
	if (mode == 0)
	{
		solve4c ();
		plot_chr ();
		plot_pot ();
	}
}

