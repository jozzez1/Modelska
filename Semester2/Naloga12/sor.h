#include <math.h>

static inline double
pow2 (double x) { return x*x; };

// get XI and its norm
void get_xi (double * xi, double * norm2, const double * psi, const double * zeta, const unsigned int N);

// the even SOR step
void step_even (double * psi, const double * w, const double * xi, unsigned int N);

// the odd SOR step
void step_odd (double * psi, const double * w, const double * xi, unsigned int N);
