#ifndef __WALKER
#define __WALKER

#include <Eigen/Dense>
#include <Eigen/LU>
#include <vector>
#include <cmath>

class Walker
{
	public:
		Walker (int);	// contructor
		~Walker ();	// destructor

		// global-local transform matrix
		Matrix2d rotMat (int);

		// calculates velocity of the flow at a certain point
		Vector2d velocity (int, Vector2d);
		
		// potention, input are global coordinates
		double potential (int, Vector2d);

		// functions to fill the matrix/vectors
		void fillu () { u = VectorXd::Ones (N-1); };
		void fillA ();

		void solve4c () { c = A.lu().solve(u) };

	private:
		int N;		// number of boundary indices
		std::vector<Vector2d> point;	// vector with all the points

		void (* init_points) (int, int);

		VectorXd c;	// vector for solving the system
		VectorXd u;	// potentials -- all of these are 1
		MatrixXd A;	// and the matrix
};

#endif
