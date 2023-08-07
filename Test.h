#include "IncludeLibraries.h"

//Test class that represents a test for the algorithm that contains a set of points and a label
class Test
{
public:
	vector<Eigen::Vector3d> points;
	bool is_wall;

	// Constructor
	Test(bool isWall_) {
		this->points = random_wall(isWall_);
		this->is_wall = isWall_;
	}

	// Sets points randomly for the test
	vector<Eigen::Vector3d> random_wall(bool is_wall);

};