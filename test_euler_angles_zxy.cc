#include "euler_angles_zxy.hpp"
#include "ltest.hpp"

#include <cmath>
#include <Eigen/Geometry>

using namespace mypilot::mymath;

Eigen::Quaterniond GoldenEulerZXYToQuaternion(const double roll,
                                              const double pitch,
                                              const double yaw) {
  return Eigen::AngleAxisd(yaw, Eigen::Vector3d::UnitZ()) *
         Eigen::AngleAxisd(pitch, Eigen::Vector3d::UnitX()) *
         Eigen::AngleAxisd(roll, Eigen::Vector3d::UnitY());
}

int main(int argc, char* argv[]) {
  TEST_START("SingleConstruct")
  {
    EulerAnglesZXYd a(1.0);
    EXPECT_DOUBLE_EQ(0.0, a.roll());
    EXPECT_DOUBLE_EQ(0.0, a.pitch());
    EXPECT_DOUBLE_EQ(1.0, a.yaw());

    auto q = a.to_quaternion();
    EXPECT_DOUBLE_EQ(cos(0.5), q.w());
    EXPECT_DOUBLE_EQ(0.0, q.x());
    EXPECT_DOUBLE_EQ(0.0, q.y());
    EXPECT_DOUBLE_EQ(sin(0.5), q.z());
  }
  TEST_END("SingleConstruct");

  TEST_START("FullConstructDouble")
  {
    EulerAnglesZXYd a(0.35, 0.24, -1.0);
    auto q_golden = GoldenEulerZXYToQuaternion(0.35, 0.24, -1.0);
    auto q = a.to_quaternion();
    EXPECT_DOUBLE_EQ(q_golden.w(), q.w());
    EXPECT_DOUBLE_EQ(q_golden.x(), q.x());
    EXPECT_DOUBLE_EQ(q_golden.y(), q.y());
    EXPECT_DOUBLE_EQ(q_golden.z(), q.z());

    EulerAnglesZXYd b(q);
    EXPECT_DOUBLE_EQ(0.35, b.roll());
    EXPECT_DOUBLE_EQ(0.24, b.pitch());
    EXPECT_DOUBLE_EQ(-1.0, b.yaw());
  }
  TEST_END("FullConstructDouble");

  TEST_START("FullConstructFloat") 
  {
    Eigen::Quaternionf q(1.0f, 2.0f, -3.0f, 4.0f);
    q.normalize();
    EulerAnglesZXYf a(q);
    auto q2 = a.to_quaternion();
    EXPECT_NEAR(q.w(), q2.w(), 5e-7);
    EXPECT_NEAR(q.x(), q2.x(), 5e-7);
    EXPECT_NEAR(q.y(), q2.y(), 5e-7);
    EXPECT_NEAR(q.z(), q2.z(), 5e-7);
  }
  TEST_END("FullConstructFloat");

  TEST_START("is_valid")
  {
    EulerAnglesZXYd a(0.35, 5.0, -1.0);
    EXPECT_TRUE(a.is_valid());

    EulerAnglesZXYd b(0.35, 2.0, -1.0);
    EXPECT_FALSE(b.is_valid());
  }
  TEST_END("is_valid");
}

