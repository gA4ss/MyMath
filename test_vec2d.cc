#include "vec2d.hpp"
#include "ltest.hpp"

#include <cmath>

using namespace mypilot::mymath;

int main(int argc, char* argv[]) {
  TEST_START("Basic");
  {
    Vec2d pt(2, 3);
    EXPECT_NEAR(pt.length(), std::sqrt(13.0), 1e-5);
    EXPECT_NEAR(pt.length_square(), 13.0, 1e-5);
    EXPECT_NEAR(pt.distance_to({0, 0}), std::sqrt(13.0), 1e-5);
    EXPECT_NEAR(pt.distance_square_to({0, 0}), 13.0, 1e-5);
    EXPECT_NEAR(pt.distance_to({0, 2}), std::sqrt(5.0), 1e-5);
    EXPECT_NEAR(pt.distance_square_to({0, 2}), 5.0, 1e-5);
    EXPECT_NEAR(pt.angle(), std::atan2(3, 2), 1e-5);
    EXPECT_NEAR(pt.cross_prod({4, 5}), -2, 1e-5);
    EXPECT_NEAR(pt.inner_prod({4, 5}), 23, 1e-5);
#ifdef MYMATH_DBG
    EXPECT_EQ(pt.str(), "vec2d ( x = 2  y = 3 )");
#endif
    pt.set_x(4);
    pt.set_y(5);
    EXPECT_NEAR(pt.length(), std::sqrt(41.0), 1e-5);
    EXPECT_NEAR(pt.length_square(), 41.0, 1e-5);
    pt.normalize();
    EXPECT_NEAR(pt.x(), 4.0 / std::sqrt(41.0), 1e-5);
    EXPECT_NEAR(pt.y(), 5.0 / std::sqrt(41.0), 1e-5);
    EXPECT_NEAR(pt.length(), 1.0, 1e-5);

    const Vec2d d = Vec2d(0.5, 1.5) + Vec2d(2.5, 3.5);
    EXPECT_NEAR(d.x(), 3.0, 1e-5);
    EXPECT_NEAR(d.y(), 5.0, 1e-5);
    const Vec2d e = Vec2d(0.5, 1.5) - Vec2d(2.5, 3.5);
    EXPECT_NEAR(e.x(), -2.0, 1e-5);
    EXPECT_NEAR(e.y(), -2.0, 1e-5);
    const Vec2d f = d / 2.0;
    EXPECT_NEAR(f.x(), 1.5, 1e-5);
    EXPECT_NEAR(f.y(), 2.5, 1e-5);
    const Vec2d g = e * (-3.0);
    EXPECT_NEAR(g.x(), 6.0, 1e-5);
    EXPECT_NEAR(g.y(), 6.0, 1e-5);

    const Vec2d unit_pt = Vec2d::create_unit_vec2d(M_PI_4);
    EXPECT_NEAR(unit_pt.x(), std::sqrt(2.0) / 2.0, 1e-5);
    EXPECT_NEAR(unit_pt.y(), std::sqrt(2.0) / 2.0, 1e-5);
    EXPECT_NEAR(unit_pt.angle(), M_PI_4, 1e-5);
  }
  TEST_END("Basic");

  TEST_START("rotate");
  {
    Vec2d pt(4, 0);
    auto p1 = pt.rotate(M_PI / 2.0);
    EXPECT_NEAR(p1.x(), 0.0, 1e-5);
    EXPECT_NEAR(p1.y(), 4.0, 1e-5);
    auto p2 = pt.rotate(M_PI);
    EXPECT_NEAR(p2.x(), -4.0, 1e-5);
    EXPECT_NEAR(p2.y(), 0.0, 1e-5);
    auto p3 = pt.rotate(-M_PI / 2.0);
    EXPECT_NEAR(p3.x(), 0.0, 1e-5);
    EXPECT_NEAR(p3.y(), -4.0, 1e-5);
    auto p4 = pt.rotate(-M_PI);
    EXPECT_NEAR(p4.x(), -4.0, 1e-5);
    EXPECT_NEAR(p4.y(), 0.0, 1e-5);
  }
  TEST_END("rotate");

  return 0;
}
