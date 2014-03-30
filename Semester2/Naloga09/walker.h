#ifndef WALKER_H
#define WALKER_H

#include <vector>
#include <cmath>
#include <iostream>
#include <string>

#include <Eigen/Dense>
#include <Eigen/QR>

#include <mgl2/mgl.h>
#include <mgl2/data.h>

typedef Eigen::Matrix2d	Matrix2d;
typedef Eigen::Vector2d	Vector2d;
typedef Eigen::MatrixXd	MatrixXd;
typedef Eigen::VectorXd	VectorXd;

class Walker
{
	public:
		// constructor
		Walker (int, int, double, double,
			int, int, double, double, double, double,
			double);

		//destructor
		~Walker ();

		// global-local transform matrix
		Matrix2d rotMat (int);

		// calculates velocity of the flow at a certain point
		Vector2d velocity (int, Vector2d);
		
		// potention, input are global coordinates
		double potential (int, Vector2d);

		// calculate capacity/force
		double effect (void);

		// functions to fill the matrix/vectors
		void init_points (int, int, double, double);
		void fillu (int);
		void fillA (int);

		// return the solution
		VectorXd solution () { return c; };

		// print out the solution
		void print_solution () { std::cout << c << std::endl; };

		// and here comes the general solver
		void solve (int);

		// plotting via MathGL
		void plot_flow ();				// flow plot for 2nd assignment
		void plot_vec (int);				// gradient plot
		void plot_pot  ();				// plot potential
		void plot_chr  ();				// plot charge density

	private:
		int N;						// number of boundary indices
		int nx;
		int ny;

		double xmax;
		double xmin;
		double ymax;
		double ymin;
		double u_inf;

		std::vector<Vector2d> point;			// vector with all the points
		std::vector<Vector2d> (* point_create) (int, void *);

		// finally, the solver itself for the 1st assignment
		void solve4c () { c = A.householderQr().solve(u); };

		// check the parallel points
		void control (int, double);

		VectorXd c;					// vector for solving the system
		VectorXd u;					// potentials -- all of these are 1
		MatrixXd A;					// and the matrix
};

#endif
