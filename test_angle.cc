#include "ltest.hpp"
#include "angle.hpp"

#include <cmath>

#ifdef USE_SIN_TABLE
#include "sin_table.h"
#endif

using namespace mypilot::mymath;

int main(int argc, char* argv[]) {

#ifdef USE_SIN_TABLE
  TEST_START("SIN_TABLE");
  {
    EXPECT_FLOAT_EQ(0.0f, SIN_TABLE[0]);
    EXPECT_FLOAT_EQ(1.0f, SIN_TABLE[16384]);
  }
  TEST_END("SIN_TABLE");
#endif

  TEST_START("Angle8");
  {
    auto a = Angle8::from_deg(90.0);
    EXPECT_DOUBLE_EQ(90.0, a.to_deg());
    EXPECT_DOUBLE_EQ(M_PI_2, a.to_rad());
    EXPECT_FLOAT_EQ(1.0f, sin(a));
    EXPECT_FLOAT_EQ(0.0f, cos(a));
  }
  TEST_END("Angle8");

  TEST_START("Angle16");
  {
    auto a = Angle16(1);
    EXPECT_DOUBLE_EQ(180.0 / 32768, a.to_deg());

    a = Angle16::from_deg(-150.0 - 360.0);
    EXPECT_NEAR(-0.5f, sin(a), 1e-4);
    EXPECT_NEAR(-0.5 * sqrt(3), cos(a), 1e-4);
  }
  TEST_END("Angle16");

  TEST_START("Angle32");
  {
    auto a = Angle32::from_rad(1.0);
    EXPECT_NEAR(180 / M_PI, a.to_deg(), 1e-7);
    EXPECT_NEAR(1.0, a.to_rad(), 1e-9);
  }
  TEST_END("Angle32");

  TEST_START("operators") {
    auto a = Angle16::from_deg(100.0);
    auto b = a;
    a += b;
    a *= 0.5;
    a = 7 * (a + b * 0.7);
    a /= 1.1;
    EXPECT_DOUBLE_EQ(-63.65478515625, a.to_deg());
  }
  TEST_END("operators");

  return 0;
}