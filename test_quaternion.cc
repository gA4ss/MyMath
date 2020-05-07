#include "quaternion.hpp"
#include "ltest.hpp"

#include <cmath>

using namespace mypilot::mymath;

int main(int argc, char* argv[]) {

  TEST_START("quaternion_to_heading");
  {
    const double v = sqrt(0.5);  // = cos(pi / 4) = sin(pi / 4)
    EXPECT_DOUBLE_EQ(0, quaternion_to_heading(v, 0.0, 0.0, -v));  // 指向东方
    EXPECT_DOUBLE_EQ(M_PI_2, quaternion_to_heading(1.0, 0.0, 0.0, 0.0));  // 指向北方
    EXPECT_DOUBLE_EQ(-M_PI_2, quaternion_to_heading(0.0, 0.0, 0.0,1.0));  // 指向南方
    EXPECT_DOUBLE_EQ(-M_PI, quaternion_to_heading(v, 0.0, 0.0, v));  // 指向西方

    Eigen::Quaternionf q(1.0, 0.0, 0.0, 0.0);
    EXPECT_FLOAT_EQ(M_PI_2, quaternion_to_heading(q));  // 指向北方

    const double headings[] = {-3.3, -2.2, -1.1, 1.2, 2.3, 3.4, 4.5, 5.6, 6.7};
    for (double heading : headings) {
      EXPECT_NEAR(normalize_angle(heading),
                  quaternion_to_heading(heading_to_quaternion<double>(heading)),
                  1e-15);
    }
  }
  TEST_END("quaternion_to_heading");

  TEST_START("quaternion_rotate") {
    Quaternion q;
    q.set_qx(0.016590540978116377);
    q.set_qy(0.012968083311103572);
    q.set_qz(-0.99256254167039326);
    q.set_qw(-0.1199007240933047);

    Eigen::Vector3d original(0.18112868882363914, 0.38614886414425986,
                            -0.15861744649897938);

    auto rotated = quaternion_rotate(q, original);
    EXPECT_NEAR(rotated[0], -0.26184808017295008, 1e-9);
    EXPECT_NEAR(rotated[1], -0.32827419468368224, 1e-9);
    EXPECT_NEAR(rotated[2], -0.17535585973456849, 1e-9);
  }
  TEST_END("quaternion_to_heading");

  TEST_START("inverse_quaternion_rotate") {
    Quaternion q;
    q.set_qx(0.016590540978116377);
    q.set_qy(0.012968083311103572);
    q.set_qz(-0.99256254167039326);
    q.set_qw(-0.1199007240933047);
    Eigen::Vector3d rotated(-0.26184808017295008, -0.32827419468368224,
                            -0.17535585973456849);
    auto original = inverse_quaternion_rotate(q, rotated);
    EXPECT_NEAR(original[0], 0.18112868882363914, 1e-9);
    EXPECT_NEAR(original[1], 0.38614886414425986, 1e-9);
    EXPECT_NEAR(original[2], -0.15861744649897938, 1e-9);
  }
  TEST_END("inverse_quaternion_rotate");
}
