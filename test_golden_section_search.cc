#include "golden_section_search.hpp"
#include "ltest.hpp"

#include <cmath>

using namespace mypilot::mymath;

double linear_func(double x) { return 2.0 * x; }
double square_func(double x) { return x * x; }
double cubic_func(double x) { return (x - 1.0) * (x - 2.0) * (x - 3.0); }
double sin_func(double x) { return std::sin(x); }

int main(int argc, char* argv[]) {
  TEST_START("golden_section_search");
  {
    double linear_argmin = golden_section_search(linear_func, 0.0, 1.0, 1e-6);
    EXPECT_NEAR(linear_argmin, 0.0, 1e-5);
    double square_argmin = golden_section_search(square_func, -1.0, 2.0, 1e-6);
    EXPECT_NEAR(square_argmin, 0.0, 1e-5);
    double cubic_argmin_1 = golden_section_search(cubic_func, 0.0, 1.5, 1e-6);
    EXPECT_NEAR(cubic_argmin_1, 0.0, 1e-5);
    double cubic_argmin_2 = golden_section_search(cubic_func, 1.0, 1.8, 1e-6);
    EXPECT_NEAR(cubic_argmin_2, 1.0, 1e-5);
    double cubic_argmin_3 = golden_section_search(cubic_func, 2.0, 3.0, 1e-6);
    EXPECT_NEAR(cubic_argmin_3, 2.0 + 1.0 / std::sqrt(3.0), 1e-5);
    double sin_argmin = golden_section_search(sin_func, 0.0, 2 * M_PI, 1e-6);
    EXPECT_NEAR(sin_argmin, 1.5 * M_PI, 1e-5);
  }
  TEST_END("golden_section_search");
}
