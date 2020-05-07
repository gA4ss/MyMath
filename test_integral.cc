#include "integral.hpp"
#include "ltest.hpp"

#include <cmath>

using namespace mypilot::mymath;

double linear_func(double x) { return 2.0 * x; }
double square_func(double x) { return x * x; }
double cubic_func(double x) { return x * x * x; }
double sin_func(double x) { return std::sin(x); }

int main(int argc, char* argv[]) {
  TEST_START("integration");
  {
    double linear_integral = integrate_by_gauss_legendre<5>(linear_func, 0.0, 1.0);
    EXPECT_NEAR(linear_integral, 1.0, 1e-5);
    double square_integral = integrate_by_gauss_legendre<5>(square_func, 0.0, 1.0);
    EXPECT_NEAR(square_integral, 1.0 / 3.0, 1e-5);
    double cubic_integral = integrate_by_gauss_legendre<5>(cubic_func, 0.0, 1.0);
    EXPECT_NEAR(cubic_integral, 1.0 / 4.0, 1e-5);
    double sin_integral = integrate_by_gauss_legendre<5>(sin_func, 0.0, 0.5 * M_PI);
    EXPECT_NEAR(sin_integral, 1.0, 1e-5);
  }
  TEST_END("integration");
}
