#include "Test.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
    //std::uniform_int_distribution<float> distribution(0, 1);

    vector<Eigen::Vector3d> random_wall;
    while (random_points > 0)
    {
        vector<double> point;
        double x = dist(gen);
        double y = dist(gen);
        double z;
        if (is_wall)
            z = gsl_ran_gaussian(rng, 0.1) + 0.5;
        else
            z = gsl_rng_uniform(rng);

        //z = (z + 3.0) / 6.0;  // Apply transformation to map it between 0 and 1
        random_wall.push_back(Eigen::Vector3d(x, y, z));
        random_points--;
    }
    //// Free GSL RNG
    gsl_rng_free(rng);

    return random_wall;
}