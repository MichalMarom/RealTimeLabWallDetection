#include "Test.h"

// Function that given if this test represent a wall or not -  Sets points randomly for that test
// input:
//       bool                   is_wall
//
// output:
//       vector<Eigen::Vector3d> random_wall (Points vector)

vector<Eigen::Vector3d> Test::random_wall(bool is_wall)
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
        z = (z + 3.0) / 6.0;  // Apply transformation to map it between 0 and 1
        random_wall.push_back(Eigen::Vector3d(x, y, z));
        random_points--;
    }

    return random_wall;
}
