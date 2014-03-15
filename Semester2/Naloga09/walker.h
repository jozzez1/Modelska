#ifndef WALKER_H
#define WALKER_H

#include <vector>
#include <cmath>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/LU>

typedef Eigen::Matrix2d	Matrix2d;
typedef Eigen::Vector2d	Vector2d;
typedef Eigen::MatrixXd	MatrixXd;
typedef Eigen::VectorXd	VectorXd;

class Walker
{
	public:
		Walker (int, int, double, double);		// contructor
		~Walker ();					// destructor

		// global-local transform matrix
		Matrix2d rotMat (int);

		// calculates velocity of the flow at a certain point
		Vector2d velocity (int, Vector2d);
		
		// potention, input are global coordinates
		double potential (int, Vector2d);

		// functions to fill the matrix/vectors
		void init_points (int, int, double, double);
		void fillu () { u = VectorXd::Ones (N-1); };
		void fillA ();

		// finally, the solver itself
		void solve4c () { c = A.lu().solve(u); };

		// return the solution
		VectorXd solution () { return c; };

		// print out the solution
		void print_solution () { std::cout << c << std::endl; };

	private:
		int N;						// number of boundary indices
		std::vector<Vector2d> point;			// vector with all the points
		std::vector<Vector2d> (* point_create) (int, void *);

		VectorXd c;					// vector for solving the system
		VectorXd u;					// potentials -- all of these are 1
		MatrixXd A;					// and the matrix
};

#endif
