#pragma once

#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <cstdlib>
#include <math.h>
#include <iomanip>
#include <Eigen/Dense>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_cdf.h>
# define M_PI	3.14159265358979323846  // Pi

using namespace std;

// WallHandle class checks if given points are a wall
class WallHandle{

public:
	bool wallDetector(vector <Eigen::Vector3d>& points);
	//bool isNormallyDistributed(const std::vector<double>& data);
	double angleBetweenPlanes(Eigen::Vector3d normal_vector);
	Eigen::Vector4d findMinimizingPlane(const vector<Eigen::Vector3d>& points);
};