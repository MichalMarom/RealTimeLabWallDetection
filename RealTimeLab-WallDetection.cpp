#include "IncludeLibraries.h"
#include "Test.h"

bool isNormallyDistributed(const std::vector<double>& data, double significance_level = 0.05) {
    if (data.empty()) {
        // Empty data vector, cannot perform the test
        return false;
    }

    // Sort the data vector in ascending order
    std::vector<double> sorted_data = data;
    std::sort(sorted_data.begin(), sorted_data.end());

    // Calculate the mean and standard deviation of the data
    double mean = gsl_stats_mean(&sorted_data[0], 1, sorted_data.size());
    double stddev = gsl_stats_sd(&sorted_data[0], 1, sorted_data.size());

    // Calculate the test statistic (D) and p-value
    double D = 0.0;
    for (size_t i = 0; i < sorted_data.size(); ++i) {
        double F_obs = gsl_cdf_ugaussian_P((sorted_data[i] - mean) / stddev);
        double F_exp = (i + 1.0) / sorted_data.size();
        double diff = std::abs(F_obs - F_exp);
        if (diff > D) {
            D = diff;
        }
    }

    // Calculate the critical value for the given significance level and sample size
    double critical_value = gsl_cdf_ugaussian_Pinv(1.0 - significance_level / 2.0) / std::sqrt(sorted_data.size());

    return D <= critical_value;
}


// Tolerance value for floating-point precision
const double tolerance = 1e-6;

vector<vector<float>> randomWall()
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

vector<vector<float>> getCollinearPoints(vector<vector<float>> points)
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

vector<float> planeEquation(float x1, float y1, float z1,
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
    /*std::cout << std::fixed;
    std::cout << std::setprecision(2);
    cout << "equation of plane is " << a << " x + " << b
        << " y + " << c << " z + " << d << " = 0.";
    */
    return equation_plane;
}

bool isPointOnPlane(vector<float>& point, vector<float>& plane)
{
    float result = (plane[0] * point[0]) + (plane[1] * point[1]) + (plane[2] * point[2]) + plane[3];
    return (abs(result) < tolerance);
}

bool planeSolution(vector <vector<float>>& wall)
{
    // Find hyperplane with 3 points from the samples and check if all the point lies in the hyperplane
    vector<vector<float>> collinear_points = getCollinearPoints(wall);
    vector<float> plane_equation = planeEquation(collinear_points[0][0], collinear_points[0][1], collinear_points[0][2],
                                                    collinear_points[1][0], collinear_points[1][1], collinear_points[1][2],
                                                    collinear_points[2][0], collinear_points[2][1], collinear_points[2][2]);
    for (auto point : wall)
    {
        if (isPointOnPlane(point, plane_equation) == false)
        {
            return false;
        }
    }
    return true;
}

///////////////////////////
// ****** T-Test ****** //
/////////////////////////
float computeMean(const vector<float> value)
{
    float sum = 0;
    for (auto e : value) { sum += e; }
    return sum / value.size();
}

float computeStdDeviation(const vector<float>& values, float mean)
{
    float varianceSum = 0;
    for (auto value : values) { varianceSum += (value - mean) * (value - mean); }
    float variance = varianceSum / values.size();
    return sqrt(variance);
}

bool TTestSol(vector <vector<float>>& wall)
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

bool zTest(const std::vector<Eigen::Vector3d>& wall)
{
    // Find the noise in Z coordinates and check if the t_score is smeller then 1 - then the noise is normally distributed
    vector<float> z_cord;
    for (const auto& point : wall)
    {
        z_cord.push_back(point[2]);
    }

    float mean = computeMean(z_cord);
    float stddev = computeStdDeviation(z_cord, mean);

    // Calculate the sample mean
    float sample_mean = computeMean(z_cord);

    // Calculate the sample standard deviation
    float sample_stddev = computeStdDeviation(z_cord, sample_mean);

    // Calculate the number of samples
    int num_samples = static_cast<int>(z_cord.size());

    // Calculate the Z-score
    double z_score = (sample_mean - mean) / (stddev / sqrt(static_cast<double>(num_samples)));

    // Define the significance level (alpha) - choose an appropriate value based on your test
    const double alpha = 0.05;

    // Check if the absolute Z-score is less than the critical value for the significance level
    double critical_value = 1.96; // For a two-tailed test at alpha = 0.05
    if (std::abs(z_score) <= critical_value)
    {
        return true; // Null hypothesis accepted, sample mean is not significantly different from population mean
    }
    else
    {
        return false; // Null hypothesis rejected, sample mean is significantly different from population mean
    }
}

bool isDataNormallyDistributed(const std::vector<Eigen::Vector3d>& wall)
{
    // Perform the Z-test on the 'wall' data
    bool is_data_normal = zTest(wall);

    // Return the result
    return is_data_normal;
}

//////////////////////////
// ****** Eigen ****** //
////////////////////////

/* Function that calculate the sum of distance to a set of points from the plane
input:
      Eigen::Vector4d&   plane
      vector<Eigen::Vector3d> points
output:       
       float error
*/
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

/* Function that find the plane that minimizes the distance to a set of points
 input:
       vector<Eigen::Vector3d> points

 output:
        Eigen::Vector4d   plane equation (ax+by+cz+d=0)
*/
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

/* Function that calculate the angle between given plane to X - Z Plane
 input:
       Eigen::Vector3d& plane (coefficients a,b,c of the plane equation)

 output:
        double    degrees_angle (in degrees)   
*/
double angleBetweenPlanes(Eigen::Vector3d normalVector)
{
    // Normal vector of the X-Z plane
    Eigen::Vector3d xzPlaneNormal = Eigen::Vector3d(0, 1, 0); 

    double dot_prod = (normalVector[0] * xzPlaneNormal[0]) + (normalVector[1] * xzPlaneNormal[1]) + (normalVector[2] * xzPlaneNormal[2]);
    double mag1 = sqrt((normalVector[0] * normalVector[0]) + (normalVector[1] * normalVector[1]) + (normalVector[2] * normalVector[2]));
    double mag2 = sqrt((xzPlaneNormal[0] * xzPlaneNormal[0]) + (xzPlaneNormal[1] * xzPlaneNormal[1]) + (xzPlaneNormal[2] * xzPlaneNormal[2]));

    double cos_angle = dot_prod / (mag1 * mag2);
    double radian_angle = acos(cos_angle);

    // Convert the angle from radians to degrees
    double degrees_angle = radian_angle * (180.0 / M_PI);
    return degrees_angle;
}

/* Function that find the XY plane
 input:
       vector<Eigen::Vector3d> points

 output:
        Eigen::Vector4d   plane equation (ax+by+cz+d=0)
*/
Eigen::Vector4d findXYPlane(const std::vector<Eigen::Vector3d>& points) {
    Eigen::Matrix<double, Eigen::Dynamic, 3> A(points.size(), 3);

    // Fill the matrix A with points
    for (size_t i = 0; i < points.size(); ++i) {
        A.row(i) = points[i].transpose();
    }

    // Compute the centroid of the points
    Eigen::Vector3d centroid = A.colwise().mean();

    // Calculate the normal vector for the XY plane (Z=0)
    Eigen::Vector3d normal(0.0, 0.0, 1.0);

    // Calculate the d value (distance from origin to the plane)
    double d = -normal.dot(centroid);

    // Return the plane as a 4D vector (a, b, c, d)
    return Eigen::Vector4d(normal.x(), normal.y(), normal.z(), d);
}

/* Function that given points - decides whether they are a wall
 input:
       vector<Eigen::Vector3d> points

 output:
        bool  is_wall
*/
bool wallDetector(vector <Eigen::Vector3d>& points)
{
    bool is_wall = false;
    plotPlane(points);
    vector<float> x_cord;
    for (auto point : points)
    {
        x_cord.push_back(point[0]);
    }
    float x_mean = computeMean(x_cord);

    vector<float> y_cord;
    for (auto point : points)
    {
        y_cord.push_back(point[1]);
    }
    float y_mean = computeMean(y_cord);

    vector<float> z_cord;
    for (auto point : points)
    {
        z_cord.push_back(point[2]);
    }

    float z_mean = computeMean(z_cord);

    float test_1 = z_mean / x_mean;
    float test_2 = z_mean / y_mean;

    test_1 = (int(floor(test_1 * 100))) / 100;
    test_2 = (int(floor(test_2 * 100))) / 100;
    int itest_1 = (int)test_1;
    int itest_2 = (int)test_2;

    test_1 -= itest_1;
    test_2 -= itest_2;

    if (test_1 == test_2)
    {
        return true;
    }
    else {
        return false;
    }

    std::vector<double> z_cord;
    for (Eigen::Vector3d point : points) {
        z_cord.push_back(point.z());
    }
    if (! isNormallyDistributed(z_cord))
    {
        is_wall = false;
        return is_wall;
    }

    // Find the plane that minimizes the distance to the points
    Eigen::Vector4d plane = findMinimizingPlane(points);

    Eigen::Vector3d plane_normal = Eigen::Vector3d(plane[0], plane[1], plane[2]);

    // Find the angle between the plane and XZ-plane
    double angle_between_plane_and_XZplane = angleBetweenPlanes(plane_normal);

    if (angle_between_plane_and_XZplane >= 88 && angle_between_plane_and_XZplane <= 92)
    {
        is_wall =  true;
    }
    else {
        is_wall = false;
    }
    return is_wall;
    
    //// Find XY-plane
    //Eigen::Vector4d XY_plane = findXYPlane(points);
    //Eigen::Vector3d XY_plane_normal = Eigen::Vector3d(XY_plane[0], XY_plane[1], XY_plane[2]);
    //// Find the angle between XY-plane and XZ-plane
    //double angle_between_XYplane_and_XZplane = angleBetweenPlanes(XY_plane_normal);

    //double isEqual = (angle_between_XYplane_and_XZplane / angle_between_plane_and_XZplane);

    //return std::abs((angle_between_XYplane_and_XZplane - angle_between_plane_and_XZplane) < 0.001);
}

int main() 
{
    std::vector<std::vector<double>> x, y, z;
    for (double i = -5; i <= 5; i += 0.25) {
        std::vector<double> x_row, y_row, z_row;
        for (double j = -5; j <= 5; j += 0.25) {
            x_row.push_back(i);
            y_row.push_back(j);
            z_row.push_back(::std::sin(::std::hypot(i, j)));
        }
        x.push_back(x_row);
        y.push_back(y_row);
        z.push_back(z_row);
    }

    //plt::plot_surface(x, y, z);
    plt::plot({ 1,3,2,4 });
    plt::show();

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
    int counter = 0;
    for (auto test : tests)
    {
        predicted_labels.push_back(wallDetector(test.points));
        counter++;
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

