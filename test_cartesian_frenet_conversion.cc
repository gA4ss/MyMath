#include "cartesian_frenet_conversion.hpp"
#include "ltest.hpp"

#include <array>
#include <cmath>

using namespace mypilot::mymath;

int main(int argc, char* argv[]) {
  TEST_START("cartesian_to_frenet_test");
  {
    double rs = 10.0;
    double rx = 0.0;
    double ry = 0.0;
    double rtheta = M_PI / 4.0;
    double rkappa = 0.1;
    double rdkappa = 0.01;
    double x = -1.0;
    double y = 1.0;
    double v = 2.0;
    double a = 0.0;
    double theta = M_PI / 3.0;
    double kappa = 0.11;

    std::array<double, 3> s_conditions;
    std::array<double, 3> d_conditions;

    CartesianFrenetConverter::cartesian_to_frenet(
        rs, rx, ry, rtheta, rkappa, rdkappa, x, y, v, a, theta, kappa,
        &s_conditions, &d_conditions);

    double x_out;
    double y_out;
    double theta_out;
    double kappa_out;
    double v_out;
    double a_out;

    CartesianFrenetConverter::frenet_to_cartesian(
        rs, rx, ry, rtheta, rkappa, rdkappa, s_conditions, d_conditions, &x_out,
        &y_out, &theta_out, &kappa_out, &v_out, &a_out);

    EXPECT_NEAR(x, x_out, 1.0e-6);
    EXPECT_NEAR(y, y_out, 1.0e-6);
    EXPECT_NEAR(theta, theta_out, 1.0e-6);
    EXPECT_NEAR(kappa, kappa_out, 1.0e-6);
    EXPECT_NEAR(v, v_out, 1.0e-6);
    EXPECT_NEAR(a, a_out, 1.0e-6);
  }
  TEST_END("cartesian_to_frenet_test");
}
