#include "IncludeLibraries.h"

class Test
{
public:
        vector<Eigen::Vector3d> points;
		bool label;

		// Constructor
		Test(bool is_wall) { 
            this->points = random_wall(is_wall);
            this->label = is_wall;
		}

        vector<Eigen::Vector3d> random_wall(bool is_wall)
        {
            int random_points = 2000;
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<float> dist(0, 1);
            std::normal_distribution<float> z_dist;

            if (is_wall) {
                std::normal_distribution<float> z_dist(0.5, 0.1);
            }
            else {
                std::uniform_real_distribution<float> z_dist(0, 1);
            }

            vector<Eigen::Vector3d> random_wall;
            while (random_points > 0)
            {
                vector<float> point;
                float x = dist(gen);
                float y = dist(gen);
                float z = z_dist(gen);
                z = (z  + 3.0) / 6.0;  // Apply transformation to map it between 0 and 1
                random_wall.push_back(Eigen::Vector3d(x, y, z));
                random_points--;
            }

            return random_wall;
        }

};



