#include "walker.h"
#include "profiles.h"

// constructor
Walker::Walker (int number, int mode, double p1, double p2,
		int Nx, int Ny, double Xmin, double Xmax, double Ymin, double Ymax,
		double uInf)
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
	u_inf	= uInf;

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
	fillu (mode);
	fillA (mode);
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
Walker::fillu (int mode)
{
	if (mode == 0)
		u = VectorXd::Ones (N-1);
	else
	{
		Vector2d uInf; uInf << u_inf, 0;
		for (int i = 0; i <= N-2; i++)
		{
			Matrix2d R = rotMat (i);
			u[i] = (R * uInf)[1];
		}
	}
}

void
Walker::fillA (int mode)
{
	if (mode == 0)
	{
		for (int i = 0; i <= N-2; i++)
		{
			for (int j = 0; j <= N-2; j++)
				A(i,j) = potential (i, (point[j] + point[j+1])/2);
		}
	}

	else
	{
		for (int i = 0; i <= N-2; i++)
		{
			for (int j = 0; j <= N-2; j++)
			{
				if (i == j)
					A (i,j) = 0.5;
				else
				{
					Matrix2d R = rotMat (j);
					Vector2d v_ij = velocity (i, (point[j] + point[j+1])/2);
					A(i,j) = (R * v_ij)[1];
				}
			}
		}
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
	gr.SetRange ('x', 0, 1);
	gr.SetRange ('y', 0, 120);
	gr.Label ('y', "\\sigma(x)");
	gr.Label ('x', "x");
	gr.Axis ();
	gr.Grid ("xy", "B;");
	gr.Plot (x, (-1)*y, "r");
	gr.Title ("Porazdelitev naboja, \\sigma(x)");
	gr.WriteEPS ("naboj.eps", "My test plot");
}

void
Walker::plot_pot (void)
{
	MatrixXd P (nx * ny, 3);

	double hx = (xmax - xmin)/(nx-1),
	       hy = (ymax - ymin)/(ny-1);

	// calculate the potential
	for (int i = 0; i <= nx-1; i++)
	{
		for (int j = 0; j <= ny-1; j++)
		{
			P (j + ny*i, 0) = hx*i + xmin;
			P (j + ny*i, 1) = hy*j + ymin;
			P (j + ny*i, 2) = 0;
			for (int k = 0; k <= N-2; k++)
				P (j + ny*i, 2) += (-1)*c(k)*potential (k, P.block(j + ny*i, 0, 1, 2).transpose());
		}
	}

	mglGraph gr;

	// create data for mgl
	mglData x (nx, ny),
		y (nx, ny),
		z (nx, ny);

	x.Set (P.block(0, 0, nx*ny, 1).data(), nx, ny);
	y.Set (P.block(0, 1, nx*ny, 1).data(), nx, ny);
	z.Set (P.block(0, 2, nx*ny, 1).data(), nx, ny);

	// and now plot it
	gr.SetSize (800, 640);
	gr.SetRanges (xmin, xmax, ymin, ymax);
	gr.SetRange ('c', z);
	gr.SetTicks ('y', 1, 5);
	gr.Colorbar ("wyqrRk>");
	gr.Axis ("x", "q");
	gr.Axis ("y", "q");
	gr.Label ('x', "x", 0);
	gr.Label ('y', "y", 0);
	gr.Label ('c', "U(x,y)", 0);
	gr.Box ("q");

	gr.Cont (x, y, z, "t");
	gr.Dens (x, y, z, "wyqrRk");
	gr.Title ("Potencial, \\phi(x,y)", "", 6);
	gr.WritePNG ("potencial.png", "Picture of the potential", true);
}

void
Walker::plot_vec ()
{
	// first we compute the velocities
	MatrixXd P (nx*nx, 4);

	double hx = (xmax - xmin)/(nx-1),
	       hy = (ymax - ymin)/(ny-1);

	// so ... the velocities ...
	for (int i = 0; i <= nx-1; i++)
	{
		for (int j = 0; j <= ny-1; j++)
		{
			P (j + ny*i, 0) = hx*i + xmin;
			P (j + ny*i, 1) = hy*j + ymin;
			P.block (j + ny*i, 2, 1, 2) << u_inf, 0;

			for (int k = 0; k <= N-2; k++)
				P.block (j + ny*i, 2, 1, 2).transpose() +=
					(-1*c(k))*velocity (k, P.block (j + ny*i, 0, 1, 2).transpose());
		}
	}

	mglData x (nx, ny),
		y (nx, ny),
		vx(nx, ny),
		vy(nx, ny),
		px(N),		// wing profile x component
		py(N);		// wing profile y component

	for (int i = 0; i <= N-1; i++)
	{
		px.a [i] = point[i][0];
		py.a [i] = point[i][1];
	}

	x.Set (P.block(0, 0, nx*ny, 1).data(), nx, ny);
	y.Set (P.block(0, 1, nx*ny, 1).data(), nx, ny);
	vx.Set(P.block(0, 2, nx*ny, 1).data(), nx, ny);
	vy.Set(P.block(0, 3, nx*ny, 1).data(), nx, ny);

	mglGraph gr;

	gr.SetSize (800, 640);
	gr.SetRanges (xmin, xmax, ymin, ymax);
	gr.SetMeshNum (35);
	gr.SetRange ('c', -1, 1);
	gr.Axis ();
	gr.Label ('x', "x");
	gr.Label ('y', "y");
	gr.Box ();

	gr.Vect (x, y, vx, vy, "fkUBbrR", "f");	// this is row-major, despite the example
//	gr.Traj (x, y, vx, vy);
//	gr.Flow (x, y, vy, vx, "vkUBbrR", "20");	// for some reason this is suddenly column-major
	for (int i = 0; i <= 19; i++)
		gr.FlowP (mglPoint (xmin, (ymax - ymin)*i/(19) + ymin), x, y, vy, vx, "ycg");
	gr.Plot (px, py, "r");
	gr.WritePNG ("vektor.png", "Velocity field");
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
	else
	{
//		std::cout << A << std::endl;
		solve4c ();
		print_solution ();
		plot_vec ();
	}
}

