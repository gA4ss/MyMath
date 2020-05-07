#include "linear_interpolation.hpp"
#include "ltest.hpp"

#include <Eigen/Dense>

using namespace mypilot::mymath;

int main(int argc, char* argv[]) {
  TEST_START("lerp 1dim");
  {
    double t0 = 0.0;
    double t1 = 1.0;
    double x0 = 2.0;
    double x1 = 4.0;

    EXPECT_NEAR(lerp(x0, t0, x1, t1, 0.4), 2.8, 1e-6);
    EXPECT_NEAR(lerp(x0, t0, x1, t1, 0.9), 3.8, 1e-6);
    EXPECT_NEAR(lerp(x0, t0, x1, t1, 1.5), 5.0, 1e-6);
    EXPECT_NEAR(lerp(x0, t0, x1, t1, -0.3), 1.4, 1e-6);
  }
  TEST_END("lerp 1dim");

  TEST_START("lerp 2dim");
  {
    double t0 = 0.0;
    double t1 = 1.0;

    Eigen::Vector2d x0(2.0, 1.0);
    Eigen::Vector2d x1(4.0, 5.0);

    Eigen::Vector2d x = lerp(x0, t0, x1, t1, 0.4);
    EXPECT_NEAR(x.x(), 2.8, 1e-6);
    EXPECT_NEAR(x.y(), 2.6, 1e-6);

    x = lerp(x0, t0, x1, t1, 1.2);
    EXPECT_NEAR(x.x(), 4.4, 1e-6);
    EXPECT_NEAR(x.y(), 5.8, 1e-6);

    x = lerp(x0, t0, x1, t1, -0.5);
    EXPECT_NEAR(x.x(), 1.0, 1e-6);
    EXPECT_NEAR(x.y(), -1.0, 1e-6);
  }
  TEST_END("lerp 2dim");

  TEST_START("slerp case1");
  {
    double t0 = 0.0;
    double t1 = 1.0;
    double a0 = -2.0;
    double a1 = 8.5;

    EXPECT_NEAR(slerp(a0, t0, a1, t1, 0.4), -2.827, 1e-3);
  }
  TEST_END("slerp case1");

  TEST_START("slerp case2");
  {
    double t0 = 0.0;
    double t1 = 1.0;
    double a0 = 3.00;
    double a1 = -3.00;

    EXPECT_NEAR(slerp(a0, t0, a1, t1, 0.5001), -3.1416, 1e-3);
  }
  TEST_END("slerp case2");
}
