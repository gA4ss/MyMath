#include "ltest.hpp"
#include "math_utils.hpp"

using namespace mypilot::mymath;

int main(int argc, char* argv[]) {

  TEST_START("cross_prod");
  {
    EXPECT_NEAR(cross_prod({0, 0}, {0, 1}, {1, 0}), -1.0, 1e-5);
    EXPECT_NEAR(cross_prod({0, 0}, {1, 0}, {0, 1}), 1.0, 1e-5);
    EXPECT_NEAR(cross_prod({0, 1}, {0, 0}, {1, 0}), 1.0, 1e-5);
    EXPECT_NEAR(cross_prod({1, 2}, {3, 4}, {5, 6}), 0.0, 1e-5);
    EXPECT_NEAR(cross_prod({1, 2}, {3, 4}, {6, 5}), -4.0, 1e-5);
    EXPECT_NEAR(cross_prod({2, 2}, {7, 5}, {3, 4}), 7.0, 1e-5);
  }
  TEST_END("cross_prod");

  TEST_START("inner_prod");
  {
    EXPECT_NEAR(inner_prod({0, 0}, {0, 1}, {1, 0}), 0.0, 1e-5);
    EXPECT_NEAR(inner_prod({0, 0}, {1, 0}, {0, 1}), 0.0, 1e-5);
    EXPECT_NEAR(inner_prod({0, 1}, {0, 0}, {1, 0}), 1.0, 1e-5);
    EXPECT_NEAR(inner_prod({1, 2}, {3, 4}, {5, 6}), 16.0, 1e-5);
    EXPECT_NEAR(inner_prod({1, 2}, {3, 4}, {6, 5}), 16.0, 1e-5);
    EXPECT_NEAR(inner_prod({2, 2}, {7, 5}, {3, 4}), 11.0, 1e-5);
    EXPECT_NEAR(inner_prod({2, 2}, {0, 0}, {3, 4}), -6.0, 1e-5);
  }
  TEST_END("inner_prod");

  TEST_START("wrap_angle");
  {
    EXPECT_NEAR(wrap_angle(-1.2), -1.2 + M_PI * 2.0, 1e-6);
    EXPECT_NEAR(wrap_angle(3.4), 3.4, 1e-6);
    EXPECT_NEAR(wrap_angle(5.6), 5.6, 1e-6);
    EXPECT_NEAR(wrap_angle(7.8), 7.8 - M_PI * 2.0, 1e-6);
    EXPECT_NEAR(wrap_angle(12.4), std::fmod(12.4, M_PI * 2.0), 1e-6);
    EXPECT_NEAR(wrap_angle(-12.4), std::fmod(-12.4, M_PI * 2.0) + M_PI * 2.0, 1e-6);
  }
  TEST_END("wrap_angle");

  TEST_START("normalize_angle");
  {
    EXPECT_DOUBLE_EQ(1.5, normalize_angle(1.5));
    EXPECT_DOUBLE_EQ(1.5 - M_PI, normalize_angle(1.5 + M_PI));
    EXPECT_DOUBLE_EQ(1.5, normalize_angle(1.5 + M_PI * 2));
    EXPECT_DOUBLE_EQ(1.5, normalize_angle(1.5 - M_PI * 2));
    EXPECT_DOUBLE_EQ(-1.5, normalize_angle(-1.5));
    EXPECT_DOUBLE_EQ(-9.0 + M_PI * 2, normalize_angle(-9.0));
    EXPECT_DOUBLE_EQ(-M_PI, normalize_angle(-M_PI));
    EXPECT_DOUBLE_EQ(-M_PI, normalize_angle(M_PI));
    EXPECT_DOUBLE_EQ(-M_PI, normalize_angle(-M_PI * 3));
    EXPECT_DOUBLE_EQ(-M_PI, normalize_angle(M_PI * 3));
    EXPECT_DOUBLE_EQ(0.0, normalize_angle(M_PI * 4));
  }
  TEST_END("normalize_angle");

  TEST_START("square");
  {
    EXPECT_DOUBLE_EQ(121.0, square(11.0));
    EXPECT_FLOAT_EQ(144.0f, square(-12.0f));
    EXPECT_EQ(169, square(-13));
    EXPECT_EQ(2147395600, square(46340));
    EXPECT_EQ(-2147479015, square(46341));  // Overflow!
  }
  TEST_END("square");

  TEST_START("sqr");
  {
    EXPECT_DOUBLE_EQ(121.0, sqr(11.0));
    EXPECT_DOUBLE_EQ(0.25, sqr(0.5));
    EXPECT_DOUBLE_EQ(169.0, sqr(-13.0));
  }
  TEST_END("sqr");

  TEST_START("sigmoid");
  {
    EXPECT_DOUBLE_EQ(0.5, sigmoid(0.0));
  }
  TEST_END("sigmoid");

  TEST_START("clamp");
  {
    EXPECT_EQ(1, clamp(1, 0, 6));
    EXPECT_EQ(6, clamp(7, 0, 6));
    EXPECT_EQ(0, clamp(-1, 0, 6));
    EXPECT_EQ(0, clamp(0, 0, 6));  // test lower bound as input
    EXPECT_EQ(6, clamp(6, 0, 6));  // test upper bound as input
  }
  TEST_END("clamp");

  TEST_START("rotate_axis");
  {
    double x, y;
    double expected_x, expected_y;
    expected_x = sqrt(2);
    expected_y = 0;
    rotate_axis(M_PI / 4, 1.0, 1.0, &x, &y);
    EXPECT_DOUBLE_EQ(expected_x, x);
    EXPECT_NEAR(expected_y, y, 1e-5);

    expected_x = 1;
    expected_y = 0;
    rotate_axis(M_PI / 2, 0.0, 1.0, &x, &y);
    EXPECT_DOUBLE_EQ(expected_x, x);
    EXPECT_NEAR(expected_y, y, 1e-5);

    expected_x = -1;
    expected_y = 0;
    rotate_axis(M_PI, 1.0, 0.0, &x, &y);
    EXPECT_DOUBLE_EQ(expected_x, x);
    EXPECT_NEAR(expected_y, y, 1e-5);
  }
  TEST_END("rotate_axis");

  return 0;
}