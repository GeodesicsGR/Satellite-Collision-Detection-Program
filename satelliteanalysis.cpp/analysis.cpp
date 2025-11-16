#include <iostream>
#include <math.h>
#include <cmath>
#include <iomanip>

using namespace std;

struct Vec3 {
    double x;
    double y;
    double z;
};

// Relative position vector, vector subtraction
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
Vec3 scalar_product(const Vec3 V1, const double t) {
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

Vec3 relative_motion(const Vec3 S1, const Vec3 S2, const Vec3 V1, const Vec3 V2, const double t) {
    return sum(diff(S1, S2), scalar_product(diff(V1, V2), 2));
}

// Angle between 2 vectors
double asd(const Vec3 a, const Vec3 b) {
    double d = dot(a, b);
    double norm_vector_a = norm_of_single_vector(a);
    double norm_vector_b = norm_of_single_vector(b);
    return acos(d/(norm_vector_a*norm_vector_b));
}

// Dot product, scalar value
double dot(const Vec3 a, const Vec3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

int main() {
    Vec3 S1 = {7000, 0, 0};
    Vec3 S2 = {7001, 1, 0};
}
