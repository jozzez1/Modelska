#ifndef PROFILES_H
#define PROFILES_H

#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include <cassert>

// just a simple line, y = 0, 0 <= x <= 1
std::vector<Eigen::Vector2d> simple (int N, void * params);

// elipsoid for the 2nd assignment
std::vector<Eigen::Vector2d> elipsoid (int N, void * params);

// now we have fish profile
std::vector<Eigen::Vector2d> fishpr (int N, void * params);

// the Zukovski wing profile
std::vector<Eigen::Vector2d> zukovski (int N, void * params);

#endif
