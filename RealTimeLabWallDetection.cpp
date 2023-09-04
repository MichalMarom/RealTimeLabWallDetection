#include "RealTimeLabWallDetection.h"
#include "Test.h"

//// **** Math Functions **** //
///*
//    Compute the mean of a vector
// Input:
//       const std::vector<double>    value
//
// Output:
//        double                      sum / double(value.size())
//*/
//float WallHandle::computeMean(const vector<float> values)
//{
//    float sum = 0;
//    for (auto value : values) 
//    { 
//        sum += value; 
//    }
// 
//    return sum / values.size();
//}
//
// /*
//Compute the Std Deviation of a vector
//Input :
//const std::vector<double>    values
//double                       mean
//
//Output :
//double                      sqrt(variance)
//* /
//float WallHandle::computeStdDeviation(const vector<float>& values, float mean)
//{
//    float variance_sum = 0;
//    for (auto value : values) 
//    { 
//        variance_sum += (value - mean) * (value - mean);
//    }
//    float variance = variance_sum / values.size();
//    return sqrt(variance);
//}

///*
//    Checks if data is normally distributed
// Input:
//       const std::vector<double>&       data
//
// Output:
//        bool                            test_statistic <= critical_value
//*/
//bool WallHandle::isNormallyDistributed(const std::vector<double>& data) 
//{
//    double significance_level = 0.05;
//    if (data.empty()) {
//        // Empty data vector, cannot perform the test
//        return false;
//    }
//
//    // Sort the data vector in ascending order
//    std::vector<double> sorted_data = data;
//    std::sort(sorted_data.begin(), sorted_data.end());
//
//    // Calculate the mean and standard deviation of the data
//    double mean = gsl_stats_mean(&sorted_data[0], 1, sorted_data.size());
//    double stddev = gsl_stats_sd(&sorted_data[0], 1, sorted_data.size());
//
//    // Calculate the test statistic and p-value
//    double test_statistic = 0.0;
//    for (size_t i = 0; i < sorted_data.size(); ++i) {
//        double F_obs = gsl_cdf_ugaussian_P((sorted_data[i] - mean) / stddev);
//        double F_exp = (i + 1.0) / sorted_data.size();
//        double diff = std::abs(F_obs - F_exp);
//        if (diff > test_statistic) {
//            test_statistic = diff;
//        }
//    }
//
//    // Calculate the critical value for the given significance level and sample size
//    double critical_value = gsl_cdf_ugaussian_Pinv(1.0 - significance_level / 2.0) / std::sqrt(sorted_data.size());
//
//    return test_statistic <= critical_value;
//}

/*
    Function that find the plane that minimizes the distance to a set of points
 Input:
       const vector<Eigen::Vector3d>&       points

 Output:
        Eigen::Vector4d   plane equation (ax+by+cz+d=0)
*/
Eigen::Vector4d WallHandle::findMinimizingPlane(const vector<Eigen::Vector3d>& points) {
    Eigen::Matrix<double, Eigen::Dynamic, 3> A(points.size(), 3);

    // Fill the matrix A with points
    for (size_t i = 0; i < points.size(); ++i) {
        A.row(i) = points[i].transpose();
    }

    // Compute the centroid of the points
    Eigen::Vector3d centroid = A.colwise().mean();

    // Subtract the centroid from each point to center them around the origin
    A.rowwise() -= centroid.transpose();

    // Compute the singular value decomposition (SVD) of A
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // The last column of V contains the normal vector of the plane
    Eigen::Vector3d normal = svd.matrixV().col(2);

    // The plane equation: ax + by + cz + d = 0
    // Compute the d value (distance from origin to the plane)
    double d = -normal.dot(centroid);

    // Return the plane as a 4D vector (a, b, c, d)
    return Eigen::Vector4d(normal.x(), normal.y(), normal.z(), d);
}

/* 
    Function that calculate the angle between given plane to X - Z Plane
 Input:
       Eigen::Vector3d&     plane (coefficients a,b,c of the plane equation)

Output:
        double              degrees_angle
*/
double WallHandle::angleBetweenPlanes(Eigen::Vector3d normal_vector)
{
    // Normal vector of the X-Z plane
    Eigen::Vector3d XZ_plane_normal = Eigen::Vector3d(0, 1, 0);

    double dot_prod = (normal_vector[0] * XZ_plane_normal[0]) + (normal_vector[1] * XZ_plane_normal[1]) + (normal_vector[2] * XZ_plane_normal[2]);
    double mag1 = sqrt((normal_vector[0] * normal_vector[0]) + (normal_vector[1] * normal_vector[1]) + (normal_vector[2] * normal_vector[2]));
    double mag2 = sqrt((XZ_plane_normal[0] * XZ_plane_normal[0]) + (XZ_plane_normal[1] * XZ_plane_normal[1]) + (XZ_plane_normal[2] * XZ_plane_normal[2]));

    double cos_angle = dot_prod / (mag1 * mag2);
    double radian_angle = acos(cos_angle);

    // Convert the angle from radians to degrees
    double degrees_angle = radian_angle * (180.0 / M_PI);
    return degrees_angle;
}

/* 
    Function that given points - decides whether they are a wall
 Input:
       vector<Eigen::Vector3d>      points

 Output:
        bool                        is_wall
*/
bool WallHandle::wallDetector(vector <Eigen::Vector3d>& points)
{
    bool is_wall = false;

    std::vector<double> z_cord;
    for (Eigen::Vector3d point : points) {
        z_cord.push_back(point.z());
    }
    /*if (!isNormallyDistributed(z_cord))
    {
        is_wall = false;
        return is_wall;
    }*/

    // Find the plane that minimizes the distance to the points
    Eigen::Vector4d plane = findMinimizingPlane(points);
    Eigen::Vector3d plane_normal = Eigen::Vector3d(plane[0], plane[1], plane[2]);

    // Find the angle between the plane and XZ-plane
    double angle_between_plane_and_XZplane = angleBetweenPlanes(plane_normal);

    if (angle_between_plane_and_XZplane >= 88 && angle_between_plane_and_XZplane <= 92)
    {
        is_wall = true;
    }
    else {
        is_wall = false;
    }
    return is_wall;
}

int main()
{
    vector<Test> tests;
    vector<bool> actual_labels;
    vector<bool> predicted_labels;
    int num_classes = 2;
    WallHandle wall;

    // Devides planes to wall and not a wall
    for (int i = 0; i < 10000; i++) {

        if (i < 5000) {
            tests.emplace_back(true);
            actual_labels.push_back(true);
        }
        else {
            tests.emplace_back(false);
            actual_labels.push_back(false);
        }
    }

    int counter = 0;

    // Check if given points represents a wall
    for (auto test : tests)
    {
        predicted_labels.push_back(wall.wallDetector(test.points));
        counter++;
    }

    vector<vector<int>> confusionMatrix(num_classes, vector<int>(num_classes, 0));

    // Update the confusion matrix based on the predicted and actual labels
    for (size_t i = 0; i < predicted_labels.size(); ++i) {
        int actual_class = actual_labels[i];
        int predicted_class = predicted_labels[i];
        confusionMatrix[actual_class][predicted_class]++;
    }

    return 0;
}
