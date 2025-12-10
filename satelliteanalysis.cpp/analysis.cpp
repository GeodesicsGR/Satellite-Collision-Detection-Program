#include <iostream>
#include <math.h>
#include <cmath>
#include <iomanip>
#include <numbers>
#include <chrono>
#include <iomanip>
#include <ctime>

using namespace std;

double M_PI = 2*acos(0.0);
using seconds_d = std::chrono::duration<double>;

struct Vec3 {
    double x;
    double y;
    double z;
};

// Vector subtraction
Vec3 diff(const Vec3 S1, const Vec3 S2) {
    Vec3 result; 
    result.x = S1.x - S2.x;
    result.y = S1.y - S2.y;
    result.z = S1.z - S2.z;
    return result;
}

// Vector addition
Vec3 sum(const Vec3 S1, const Vec3 S2) {
    Vec3 result;
    result.x = S1.x + S2.x;
    result.y = S1.y + S2.y;
    result.z = S1.z + S2.z;
    return result;
}

// Scalar multiplication    
Vec3 scalar_product(const Vec3 V1, const double t)      {
    Vec3 result;
    result.x = V1.x*t;
    result.y = V1.y*t;
    result.z = V1.z*t;
    return result;
}

// Norm of 2 vectors, magnitude
double norm_of_vector(Vec3 S1, Vec3 S2) {   
    Vec3 d = diff(S1, S2);
    return sqrt(pow(d.x, 2) + pow(d.y, 2) + pow(d.z, 2));
}

// Norm of a single vector, magnitude
double norm_of_single_vector(const Vec3 S1) {
    return sqrt(S1.x*S1.x + S1.y*S1.y + S1.z*S1.z);
}

// Relative motion, linear approximation
Vec3 relative_motion(const Vec3 S1, const Vec3 S2, const Vec3 V1, const Vec3 V2, const double t) {
    return sum(diff(S1, S2), scalar_product(diff(V1, V2), t));
}

// Dot product, scalar value
double dot(const Vec3 a, const Vec3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}   

// Angle between 2 vectors
double asd(const Vec3 a, const Vec3 b) {
    double d = dot(a, b);
    double norm_vector_a = norm_of_single_vector(a);
    double norm_vector_b = norm_of_single_vector(b);
    return acos(d/(norm_vector_a*norm_vector_b));
}

// Cross product
Vec3 cross_product(const Vec3 S1, const Vec3 V1) {
    Vec3 result;
    result.x = S1.y*V1.z - S1.z*V1  .y;
    result.y = S1.z*V1.x - S1.x*V1.z;
    result.z = S1.x*V1.y - S1.y*V1.x;
    return result;
}

// Unit vector along the z-axis
Vec3 k = {0, 0, 1};

// Angular momentum vector
Vec3 amv(const Vec3 S1, const Vec3 V1) {
    return cross_product(S1, V1);
}

// Node vector
Vec3 node_vector(const Vec3 S1, const Vec3 V1) {
    Vec3 amv = cross_product(S1, V1);
    return cross_product(amv, k);
}

// Orbital inclination
double orbital_inclination(const Vec3 S1, const Vec3 V1) {
    Vec3 h = amv(S1, V1);
    return asd(h, k);
}

// Right Accession of the Ascending Node (RAAN)
double raan(const Vec3 S1, const Vec3 V1, const Vec3 h, const Vec3 n) {
    double result;
    double mag = norm_of_single_vector(n);
    const double epsilon = 1e-9;
    if (mag < epsilon)
    return NAN;
    if (n.y > 0) {
        double v = (n.x/mag);
        if (v > 1.0) v = 1.0;
        if (v < -1.0) v = 1.0;

        double angle = acos(v);

        if (n.y >= 0.0) {
            return angle;
        } else {
            return 2.0*M_PI - angle;
        }
    }
}

const double EARTH_SEMI_MAJOR_AXIS_METERS = 149598023000.0;

double mean_motion() {
    return sqrt(GM/EARTH_SEMI_MAJOR_AXIS_METERS);
}

// Mean anomaly at time t
double calculate_mat(const double t, double const M0, const double T0) {
    
}

// Specific angular momentum, node vector, eccentricity vector, specific orbital energy, semi-major axis, eccentricity, inclination, RAAN, argument of perigee, true anomaly, eccentric anomaly,
// mean anomaly at epoch, mean motion (rad/s)

// Gravitational constant
constexpr double G = 6.67430e-11;

// Earth's standard gravitational parameter
constexpr double GM = 3.986004418e14;

int main() {
    Vec3 S1 = {7000, 0, 0};
    Vec3 S2 = {7001, 1, 0};
    Vec3 V1 = {0, 7.5, 0};
    Vec3 V2 = {0, 7.4, 0};

    double t = 10; // seconds
    Vec3 relative_position = relative_motion(S1, S2, V1, V2, t);
    cout << "Relative position: "
         << relative_position.x << ","
         << relative_position.y << ","
         << relative_position.z << endl;
}
