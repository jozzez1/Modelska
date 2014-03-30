#include "global.h"

unsigned int
dotplus (unsigned int m, unsigned int n)
{
    return (4*(m+n)/3) & 3;
}

double
pow2 (double a)
{
    return a*a;
}

double
dist2 (point3 a, point3 b)
{
    return sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y) * (a.z - b.z)*(a.z - b.z));
}

