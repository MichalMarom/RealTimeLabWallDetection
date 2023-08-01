#include "IncludeLibraries.h"
#include "Test.h"

// Tolerance value for floating-point precision
const double tolerance = 1e-6;

vector<vector<float>> random_wall()
{
    float random_points = 400;
    int min_cord = 1;
    float max_cord = 20;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(min_cord, max_cord);
    std::uniform_real_distribution<> z_dist(0, 1);

    vector<vector<float>> random_wall;

    while (random_points > 0)
    {
        vector<float> point;
        float x = dist(gen);
        float y = dist(gen);
        float z = z_dist(gen);
        point.push_back(x);
        point.push_back(y);
        point.push_back(z);
        random_wall.push_back(point);
        random_points--;
    }

    return random_wall;
}

///////////////////////////////
// ****** Find Plane ****** //
/////////////////////////////
bool collinear(float x1, float y1, float x2, float y2, float x3, float y3) 
{
    return (y1 - y2) * (x1 - x3) == (y1 - y3) * (x1 - x2);
}

vector<vector<float>> get_collinear_points(vector<vector<float>> points)
{
    bool is_collinear = false;
    int a_index = 0;
    int b_index = 0;
    int c_index = 0;
    srand(time(0));
    vector<vector<float>> collinear_points;

    while (is_collinear == false)
    {
        a_index = 0 + (rand() % points.size());
        b_index = 0 + (rand() % points.size());
        c_index = 0 + (rand() % points.size());
        is_collinear = collinear(points[a_index][0],
            points[a_index][1],
            points[b_index][0],
            points[b_index][1],
            points[c_index][0],
            points[c_index][1]);
    }
    collinear_points.push_back(points[a_index]);
    collinear_points.push_back(points[b_index]);
    collinear_points.push_back(points[c_index]);

    return collinear_points;
}

vector<float> equation_plane(float x1, float y1, float z1,
    float x2, float y2, float z2,
    float x3, float y3, float z3)
{
    vector<float> equation_plane;
    float a1 = x2 - x1;
    float b1 = y2 - y1;
    float c1 = z2 - z1;
    float a2 = x3 - x1;
    float b2 = y3 - y1;
    float c2 = z3 - z1;
    float a = b1 * c2 - b2 * c1;
    float b = a2 * c1 - a1 * c2;
    float c = a1 * b2 - b1 * a2;
    float d = (-a * x1 - b * y1 - c * z1);

    equation_plane.push_back(a);
    equation_plane.push_back(b);
    equation_plane.push_back(c);
    equation_plane.push_back(d);
    //std::cout << std::fixed;
    //std::cout << std::setprecision(2);
    //cout << "equation of plane is " << a << " x + " << b
    //    << " y + " << c << " z + " << d << " = 0.";

    return equation_plane;


}


bool is_point_on_plane(vector<float>& point, vector<float>& plane)
{
    float result = (plane[0] * point[0]) + (plane[1] * point[1]) + (plane[2] * point[2]) + plane[3];
    return (abs(result) < tolerance);
}

bool Sol_with_plane(vector <vector<float>>& wall)
{
    //Find hyperplane with 3 points from the samples and check if all the point lies in the hyperplane
    vector<vector<float>> collinear_points = get_collinear_points(wall);
    vector<float> plane_equation = equation_plane(collinear_points[0][0], collinear_points[0][1], collinear_points[0][2],
                                                    collinear_points[1][0], collinear_points[1][1], collinear_points[1][2],
                                                    collinear_points[2][0], collinear_points[2][1], collinear_points[2][2]);
    for (auto point : wall)
    {
        if (is_point_on_plane(point, plane_equation) == false)
        {
            return false;
        }
    }
    return true;
}

///////////////////////////
// ****** T-Test ****** //
/////////////////////////
float computeMean(const vector<float> v)
{
    float sum = 0;
    for (auto e : v) { sum += e; }
    return sum / v.size();
}

float computeStdDeviation(const vector<float>& values, float mean)
{
    float varianceSum = 0;
    for (auto value : values) { varianceSum += (value - mean) * (value - mean); }
    float variance = varianceSum / values.size();
    return sqrt(variance);
}

bool Sol_with_t_test(vector <vector<float>>& wall)
{
    // Find the noise in Z coordinates and check if the t_score is smeller then 1 - then the noise is normally distributed
    vector<float> z_cord;
    for (auto point : wall)
    {
        z_cord.push_back(point[2]);
    }
    float mean = computeMean(z_cord);
    float std = computeStdDeviation(z_cord, mean);

    int num_samples = 30;
    vector<float> z_cord_samples(num_samples);
    copy(z_cord.begin(), z_cord.begin() + num_samples, z_cord_samples.begin());
    float samples_mean = computeMean(z_cord_samples);

    float t_score = abs((sqrt(num_samples) * (mean - samples_mean)) / std);

    if (t_score > 1)
    {
        return false;
    }

    return true;
}


//////////////////////////
// ****** Eigen ****** //
////////////////////////

// Function that calculate the sum of distance to a set of points from the plane
// input:
//       Eigen::Vector4d&        plane
//       vector<Eigen::Vector3d> points
//
// output:
//        float                   error

float findPlaneError(const Eigen::Vector4d& plane, const vector<Eigen::Vector3d>& points)
{
    float error = 0;
    float a = plane[0];
    float b = plane[1];
    float c = plane[2];
    float d = plane[3];

    for (auto point : points) {
        float x = point[0];
        float y = point[1];
        float z = point[2];
        float numerator = std::abs(a * x + b * y + c * z + d);
        float denominator = std::sqrt(a * a + b * b + c * c);
        error += (numerator / denominator);

    }

    return error;
}

// Function that find the plane that minimizes the distance to a set of points
// input:
//       vector<Eigen::Vector3d> points
//
// output:
//        Eigen::Vector4d        plane equation (ax+by+cz+d=0)

Eigen::Vector4d findMinimizingPlane(const vector<Eigen::Vector3d>& points) {
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

// Function that calculate the angle between given plane to X-Z Plane
// input:
//       Eigen::Vector3d& plane (coefficients a,b,c of the plane equation)
//
// output:
//        double          angleDeg (in degrees)   

double angleBetweenPlanes(Eigen::Vector3d normalVector)
{
    // Normal vector of the X-Z plane
    Eigen::Vector3d xzPlaneNormal = Eigen::Vector3d(0, 1, 0); 

    double dotProd = (normalVector[0] * xzPlaneNormal[0]) + (normalVector[1] * xzPlaneNormal[1]) + (normalVector[2] * xzPlaneNormal[2]);
    double mag1 = sqrt((normalVector[0] * normalVector[0]) + (normalVector[1] * normalVector[1]) + (normalVector[2] * normalVector[2]));
    double mag2 = sqrt((xzPlaneNormal[0] * xzPlaneNormal[0]) + (xzPlaneNormal[1] * xzPlaneNormal[1]) + (xzPlaneNormal[2] * xzPlaneNormal[2]));

    double cosAngle = dotProd / (mag1 * mag2);
    double angleRad = acos(cosAngle);

    // Convert the angle from radians to degrees
    double angleDeg = angleRad * (180.0 / M_PI);
    return angleDeg;
}

// Function that given points - decides whether they are a wall
// input:
//       vector<Eigen::Vector3d> points
//
// output:
//        bool                   is_wall

bool wall_detector(vector <Eigen::Vector3d>& points)
{
    bool is_wall = false;

    // Find the plane that minimizes the distance to the points
    Eigen::Vector4d plane = findMinimizingPlane(points);
    Eigen::Vector3d plane_normal = Eigen::Vector3d(plane[0], plane[1], plane[2]);

    // Find the angle between the planeand XZ-plane
    double angle_between_plane_and_XZplane = angleBetweenPlanes(plane_normal);

    // Check if it is wall ??

    return true;

}

int main() 
{
    vector<Test> tests;
    vector<bool> actual_labels;
    vector<bool> predicted_labels;
    int numClasses = 2;

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


    for (auto test : tests)
    {
        predicted_labels.push_back(wall_detector(test.points));
    }

    vector<vector<int>> confusionMatrix(numClasses, vector<int>(numClasses, 0));

    // Update the confusion matrix based on the predicted and actual labels
    for (size_t i = 0; i < predicted_labels.size(); ++i) {
        int actual_class = actual_labels[i];
        int predicted_class = predicted_labels[i];
        confusionMatrix[actual_class][predicted_class]++;
    }

    return 0;
   
}

