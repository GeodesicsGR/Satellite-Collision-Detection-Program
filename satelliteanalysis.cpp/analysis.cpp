#include <iostream> // std;cout, std::end1
#include <cmath>    // sqrt, fabs
#include <iomanip>  // std::setprecision

using namespace std;

struct Vec3 {
    double x;
    double y;
    double z;
};

Vec3 subtract(const Vec3& S1, const Vec3& S2) {
    Vec3 result;
    result.x = S1.x - S2.x;
    result.y = S1.y - S2.y;
    result.z = S1.z - S2.z;
    return result;
};

int main() {
    Vec3 S1 = {7000, 0, 0};
    Vec3 S2 = {7001, 1, 0};
};