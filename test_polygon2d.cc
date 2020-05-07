#include "polygon2d.hpp"
#include "ltest.hpp"

#include <algorithm>

#ifdef MYMATH_DBG
#include <string>
#include <sstream>
#endif

#include "vec2d.hpp"
#include "box2d.hpp"
#include "line_segment2d.hpp"
#include "math_utils.hpp"

using namespace mypilot::mymath;

bool project_by_x_slow(const std::vector<Vec2d> &points, double x,
                    double *const min_y, double *const max_y) {
  *min_y = std::numeric_limits<double>::infinity();
  *max_y = -std::numeric_limits<double>::infinity();
  for (const Vec2d &p1 : points) {
    if (p1.x() < x) {
      for (const Vec2d &p2 : points) {
        if (p2.x() > x) {
          const double y = ((p2.x() - x) * p1.y() + (x - p1.x()) * p2.y()) /
                           (p2.x() - p1.x());
          *min_y = std::min(*min_y, y);
          *max_y = std::max(*max_y, y);
        }
      }
    }
  }
  return *min_y <= *max_y;
}

int main(int argc, char* argv[]) {

  TEST_START("is_point_in");
  {
    const Polygon2d poly1(Box2d::create_aabox({0, 0}, {1, 1}));
#ifdef MYMATH_DBG 
    EXPECT_EQ(poly1.str(),
              "polygon2d (  num_points = 4  points = ( ( 1, 0 ),"
              "( 1, 1 ),( 0, 1 ),( 0, 0 ) )  "
              "convex  area = 1 )");
#endif
    EXPECT_TRUE(poly1.is_point_in({0.5, 0.5}));
    EXPECT_TRUE(poly1.is_point_in({0.2, 0.2}));
    EXPECT_TRUE(poly1.is_point_in({0.2, 0.8}));
    EXPECT_TRUE(poly1.is_point_in({0.8, 0.2}));
    EXPECT_TRUE(poly1.is_point_in({0.8, 0.8}));

    EXPECT_TRUE(poly1.is_point_on_boundary({0.0, 0.0}));
    EXPECT_TRUE(poly1.is_point_in({0.0, 0.0}));
    EXPECT_TRUE(poly1.is_point_on_boundary({0.0, 0.5}));
    EXPECT_TRUE(poly1.is_point_in({0.0, 0.5}));
    EXPECT_TRUE(poly1.is_point_on_boundary({1.0, 0.5}));
    EXPECT_TRUE(poly1.is_point_in({1.0, 0.5}));

    EXPECT_FALSE(poly1.is_point_in({-0.2, 0.5}));
    EXPECT_FALSE(poly1.is_point_in({1.2, 0.5}));
    EXPECT_FALSE(poly1.is_point_in({0.5, -0.2}));
    EXPECT_FALSE(poly1.is_point_in({0.5, 1.2}));

    EXPECT_FALSE(poly1.is_point_in({0, -0.1}));
    EXPECT_FALSE(poly1.is_point_in({-0.1, 0}));
    EXPECT_FALSE(poly1.is_point_in({1.0, -0.1}));
    EXPECT_FALSE(poly1.is_point_in({-0.1, 1.0}));
    EXPECT_FALSE(poly1.is_point_in({0, 1.1}));
    EXPECT_FALSE(poly1.is_point_in({1.1, 0}));
    EXPECT_FALSE(poly1.is_point_in({1.0, 1.1}));
    EXPECT_FALSE(poly1.is_point_in({1.1, 1.0}));

    const Polygon2d poly2({{0, 1}, {1, 0}, {0, -1}, {-1, 0}});
    EXPECT_TRUE(poly2.is_point_in({0, 0}));
    EXPECT_TRUE(poly2.is_point_in({0, 0.9}));
    EXPECT_TRUE(poly2.is_point_in({0.9, 0}));
    EXPECT_TRUE(poly2.is_point_in({0, -0.9}));
    EXPECT_TRUE(poly2.is_point_in({-0.9, 0}));

    EXPECT_FALSE(poly2.is_point_in({0, 1.1}));
    EXPECT_FALSE(poly2.is_point_in({1.1, 0}));
    EXPECT_FALSE(poly2.is_point_in({0, -1.1}));
    EXPECT_FALSE(poly2.is_point_in({-1.1, 0}));

    const Polygon2d poly3({{4, 4}, {5, 6}, {6, 6}});
    EXPECT_FALSE(poly3.is_point_in({5, 4.5}));
    EXPECT_FALSE(poly3.is_point_on_boundary({5, 4.5}));
    EXPECT_TRUE(poly3.is_point_in({5, 5}));
    EXPECT_TRUE(poly3.is_point_on_boundary({5, 5}));
    EXPECT_TRUE(poly3.is_point_in({5, 5.5}));
    EXPECT_FALSE(poly3.is_point_on_boundary({5, 5.5}));
    EXPECT_TRUE(poly3.is_point_in({5, 6}));
    EXPECT_TRUE(poly3.is_point_on_boundary({5, 6}));
    EXPECT_FALSE(poly3.is_point_in({5, 6.5}));
    EXPECT_FALSE(poly3.is_point_on_boundary({5, 6.5}));

    // Concave polygons.
    const Polygon2d poly4({{0, 0}, {2, 0}, {2, 2}, {1, 1}, {0, 2}});
    EXPECT_TRUE(poly4.is_point_in({0.5, 1.5}));
    EXPECT_TRUE(poly4.is_point_on_boundary({0.5, 1.5}));
    EXPECT_FALSE(poly4.is_point_in({1.0, 1.5}));
    EXPECT_TRUE(poly4.is_point_in({1.5, 1.5}));
    EXPECT_TRUE(poly4.is_point_on_boundary({1.5, 1.5}));

    const Polygon2d poly5(
        {{0, 0}, {4, 0}, {4, 2}, {3, 2}, {2, 1}, {1, 2}, {0, 2}});
    EXPECT_FALSE(poly5.is_point_in({-0.5, 2.0}));
    EXPECT_TRUE(poly5.is_point_in({0.5, 2.0}));
    EXPECT_FALSE(poly5.is_point_in({1.5, 2.0}));
    EXPECT_FALSE(poly5.is_point_in({2.0, 2.0}));
    EXPECT_FALSE(poly5.is_point_in({2.5, 2.0}));
    EXPECT_TRUE(poly5.is_point_in({3.5, 2.0}));
    EXPECT_FALSE(poly5.is_point_in({4.5, 2.0}));
  }
  TEST_END("is_point_in");

  TEST_START("distance_to_point");
  {
    const Box2d box1(Box2d::create_aabox({0, 0}, {1, 1}));
    const Polygon2d poly1(box1);
    EXPECT_NEAR(poly1.distance_to({0.5, 0.5}), 0.0, 1e-5);

    EXPECT_NEAR(poly1.distance_to({-0.2, 0.5}), 0.2, 1e-5);
    EXPECT_NEAR(poly1.distance_to({1.2, 0.5}), 0.2, 1e-5);
    EXPECT_NEAR(poly1.distance_to({0.5, -0.2}), 0.2, 1e-5);
    EXPECT_NEAR(poly1.distance_to({0.5, 1.2}), 0.2, 1e-5);

    EXPECT_NEAR(poly1.distance_to({0, -0.1}), 0.1, 1e-5);
    EXPECT_NEAR(poly1.distance_to({-0.1, 0}), 0.1, 1e-5);
    EXPECT_NEAR(poly1.distance_to({1.0, -0.1}), 0.1, 1e-5);
    EXPECT_NEAR(poly1.distance_to({-0.1, 1.0}), 0.1, 1e-5);
    EXPECT_NEAR(poly1.distance_to({0, 1.1}), 0.1, 1e-5);
    EXPECT_NEAR(poly1.distance_to({1.1, 0}), 0.1, 1e-5);
    EXPECT_NEAR(poly1.distance_to({1.0, 1.1}), 0.1, 1e-5);
    EXPECT_NEAR(poly1.distance_to({1.1, 1.0}), 0.1, 1e-5);

    EXPECT_NEAR(poly1.distance_to({-0.1, -0.1}), 0.1 * sqrt(2.0), 1e-5);
    EXPECT_NEAR(poly1.distance_to({-0.1, 1.1}), 0.1 * sqrt(2.0), 1e-5);
    EXPECT_NEAR(poly1.distance_to({1.1, -0.1}), 0.1 * sqrt(2.0), 1e-5);
    EXPECT_NEAR(poly1.distance_to({1.1, 1.1}), 0.1 * sqrt(2.0), 1e-5);

    for (int iter = 0; iter < 10000; ++iter) {
      const double x = random_double(-10, 10);
      const double y = random_double(-10, 10);
      EXPECT_NEAR(poly1.distance_to({x, y}), box1.distance_to({x, y}), 1e-5);
    }
    for (int iter = 0; iter < 100; ++iter) {
      const double center_x = random_double(-10, 10);
      const double center_y = random_double(-10, 10);
      const double heading = random_double(0, M_PI * 2.0);
      const double length = random_double(1, 5);
      const double width = random_double(1, 5);
      const Box2d box({center_x, center_y}, heading, length, width);
      const Polygon2d polygon(box);
      for (int iter2 = 0; iter2 < 100; ++iter2) {
        const double x = random_double(-20, 20);
        const double y = random_double(-20, 20);
        EXPECT_NEAR(polygon.distance_to({x, y}), box.distance_to({x, y}), 1e-5);
      }
    }

    const Polygon2d poly2({{0, 1}, {1, 0}, {0, -1}, {-1, 0}});
    EXPECT_NEAR(poly2.distance_to({0, 0}), 0.0, 1e-5);

    EXPECT_NEAR(poly2.distance_to({0, 1.1}), 0.1, 1e-5);
    EXPECT_NEAR(poly2.distance_to({1.1, 0}), 0.1, 1e-5);
    EXPECT_NEAR(poly2.distance_to({0, -1.1}), 0.1, 1e-5);
    EXPECT_NEAR(poly2.distance_to({-1.1, 0}), 0.1, 1e-5);

    EXPECT_NEAR(poly2.distance_to({0.5, 0.5}), 0.0, 1e-5);
    EXPECT_NEAR(poly2.distance_to({0.6, 0.6}), 0.1 * sqrt(2.0), 1e-5);

    const Polygon2d poly3(
        {{0, 0}, {4, 0}, {4, 2}, {3, 2}, {2, 1}, {1, 2}, {0, 2}});
    EXPECT_NEAR(poly3.distance_to({-0.5, 2.0}), 0.5, 1e-5);
    EXPECT_NEAR(poly3.distance_to({0.5, 2.0}), 0.0, 1e-5);
    EXPECT_NEAR(poly3.distance_to({1.5, 2.0}), 0.5 / sqrt(2.0), 1e-5);
    EXPECT_NEAR(poly3.distance_to({2.0, 2.0}), 1.0 / sqrt(2.0), 1e-5);
    EXPECT_NEAR(poly3.distance_to({2.5, 2.0}), 0.5 / sqrt(2.0), 1e-5);
    EXPECT_NEAR(poly3.distance_to({3.5, 2.0}), 0.0, 1e-5);
    EXPECT_NEAR(poly3.distance_to({4.5, 2.0}), 0.5, 1e-5);
    EXPECT_NEAR(poly3.distance_to({-0.5, 1.0}), 0.5, 1e-5);
    EXPECT_NEAR(poly3.distance_to({2.0, 1.0}), 0.0, 1e-5);
    EXPECT_NEAR(poly3.distance_to({4.5, 1.0}), 0.5, 1e-5);
  }
  TEST_END("distance_to_point");

  TEST_START("distance_to_lineSegment");
  {
    const Polygon2d poly1(Box2d::create_aabox({0, 0}, {1, 1}));
    EXPECT_NEAR(poly1.distance_to({{0.5, 0.5}, {1.0, 1.0}}), 0.0, 1e-5);
    EXPECT_NEAR(poly1.distance_to({{-0.2, 0.5}, {1.2, 0.5}}), 0.0, 1e-5);
    EXPECT_NEAR(poly1.distance_to({{-2.0, -2.0}, {2.0, 2.0}}), 0.0, 1e-5);

    EXPECT_NEAR(poly1.distance_to({{-0.2, 0.5}, {-0.2, 0.8}}), 0.2, 1e-5);
    EXPECT_NEAR(poly1.distance_to({{1.2, 0.5}, {1.2, 0.3}}), 0.2, 1e-5);
    EXPECT_NEAR(poly1.distance_to({{0.5, -0.2}, {0.8, -0.2}}), 0.2, 1e-5);
    EXPECT_NEAR(poly1.distance_to({{0.5, 1.2}, {0.3, 1.2}}), 0.2, 1e-5);

    EXPECT_NEAR(poly1.distance_to({{-0.3, 0.5}, {-0.2, 0.8}}), 0.2, 1e-5);
    EXPECT_NEAR(poly1.distance_to({{1.2, 0.5}, {1.3, 0.3}}), 0.2, 1e-5);
    EXPECT_NEAR(poly1.distance_to({{0.5, -0.3}, {0.8, -0.2}}), 0.2, 1e-5);
    EXPECT_NEAR(poly1.distance_to({{0.5, 1.2}, {0.3, 1.3}}), 0.2, 1e-5);

    EXPECT_NEAR(poly1.distance_to({{0, -0.1}, {-0.1, 0}}), 0.1 / sqrt(2.0), 1e-5);
    EXPECT_NEAR(poly1.distance_to({{1.0, 1.1}, {1.1, 1.0}}), 0.1 / sqrt(2.0),
                1e-5);

    EXPECT_NEAR(poly1.distance_to({{-10.0, 0.5}, {2.0, 0.5}}), 0.0, 1e-5);
    EXPECT_NEAR(poly1.distance_to({{-1.0, 0.5}, {10.0, 0.5}}), 0.0, 1e-5);

    EXPECT_NEAR(poly1.distance_to({{-1.0, 2.0}, {-1.0, 2.0}}), sqrt(2.0), 1e-5);
    EXPECT_NEAR(poly1.distance_to({{0.5, 0.5}, {0.5, 0.5}}), 0.0, 1e-5);

    const Polygon2d poly2({{0, 1}, {1, 0}, {0, -1}, {-1, 0}});
    EXPECT_NEAR(poly2.distance_to({{-2, 0}, {2, 0}}), 0.0, 1e-5);
    EXPECT_NEAR(poly2.distance_to({{0, -2}, {0, 2}}), 0.0, 1e-5);

    EXPECT_NEAR(poly2.distance_to({{0, 1.1}, {1.1, 0}}), 0.1 / sqrt(2.0), 1e-5);
    EXPECT_NEAR(poly2.distance_to({{0, 1.1}, {-1.1, 0}}), 0.1 / sqrt(2.0), 1e-5);
    EXPECT_NEAR(poly2.distance_to({{0, -1.1}, {1.1, 0}}), 0.1 / sqrt(2.0), 1e-5);
    EXPECT_NEAR(poly2.distance_to({{0, -1.1}, {-1.1, 0}}), 0.1 / sqrt(2.0), 1e-5);
    EXPECT_NEAR(poly2.distance_to({{0.6, 0.6}, {0.7, 0.7}}), 0.1 * sqrt(2.0),
                1e-5);
    EXPECT_NEAR(poly2.distance_to({{-0.6, -0.6}, {-0.7, -0.7}}), 0.1 * sqrt(2.0),
                1e-5);
    EXPECT_NEAR(poly2.distance_to({{-0.6, -0.6}, {0.7, 0.7}}), 0.0, 1e-5);

    const Polygon2d poly3({{0, 0}, {2, 0}, {2, 2}, {1, 1}, {0, 2}});
    EXPECT_NEAR(poly3.distance_to({{-2, 0}, {2, 0}}), 0.0, 1e-5);
    EXPECT_NEAR(poly3.distance_to({{0.7, 2.0}, {1.2, 2.0}}), 0.7 / sqrt(2.0),
                1e-5);
    EXPECT_NEAR(poly3.distance_to({{0.7, 2.0}, {1.4, 2.0}}), 0.6 / sqrt(2.0),
                1e-5);
    EXPECT_NEAR(poly3.distance_to({{0.7, 1.5}, {1.6, 1.5}}), 0.0, 1e-5);
  }
  TEST_END("distance_to_lineSegment");

  TEST_START("distance_to_polygon");
  {
    const Polygon2d poly1(Box2d::create_aabox({0, 0}, {1, 1}));
    const Polygon2d poly2({{0, 1}, {1, 0}, {0, -1}, {-1, 0}});
    const Polygon2d poly3(Box2d::create_aabox({2, 2}, {3, 3}));
    const Polygon2d poly4(Box2d::create_aabox({-10, -10}, {10, 10}));

    EXPECT_NEAR(poly1.distance_to(poly2), 0.0, 1e-5);
    EXPECT_NEAR(poly1.distance_to(poly3), sqrt(2.0), 1e-5);
    EXPECT_NEAR(poly2.distance_to(poly3), 1.5 * sqrt(2.0), 1e-5);
    EXPECT_NEAR(poly4.distance_to(poly1), 0.0, 1e-5);
    EXPECT_NEAR(poly4.distance_to(poly2), 0.0, 1e-5);
    EXPECT_NEAR(poly4.distance_to(poly3), 0.0, 1e-5);
  }
  TEST_END("distance_to_polygon");

  TEST_START("contain polygon");
  {
    const Polygon2d poly1(Box2d::create_aabox({0, 0}, {3, 3}));
    const Polygon2d poly2(Box2d::create_aabox({1, 1}, {2, 2}));
    const Polygon2d poly3(Box2d::create_aabox({1.5, 1.5}, {4, 4}));
    const Polygon2d poly4(Box2d::create_aabox({-10, -10}, {10, 10}));
    EXPECT_TRUE(poly1.contains(poly2));
    EXPECT_FALSE(poly2.contains(poly1));

    EXPECT_FALSE(poly1.contains(poly3));
    EXPECT_FALSE(poly2.contains(poly3));
    EXPECT_FALSE(poly3.contains(poly1));
    EXPECT_FALSE(poly3.contains(poly2));

    EXPECT_TRUE(poly4.contains(poly1));
    EXPECT_TRUE(poly4.contains(poly2));
    EXPECT_TRUE(poly4.contains(poly3));

    const Polygon2d poly5(
        {{0, 0}, {4, 0}, {4, 2}, {3, 2}, {2, 1}, {1, 2}, {0, 2}});
    const Polygon2d poly6({{0, 1}, {4, 1}, {4, 2}, {0, 2}});
    const Polygon2d poly7({{0, 1}, {1, 1}, {1, 2}, {0, 2}});
    const Polygon2d poly8({{3, 1}, {4, 1}, {4, 2}, {3, 2}});
    const Polygon2d poly9({{0, 0}, {4, 0}, {4, 1}, {0, 1}});
    EXPECT_FALSE(poly5.contains(poly6));
    EXPECT_TRUE(poly5.contains(poly7));
    EXPECT_TRUE(poly5.contains(poly8));
    EXPECT_TRUE(poly5.contains(poly9));
  }
  TEST_END("contain polygon");

  TEST_START("convex hull");
  {
    Polygon2d polygon;
    EXPECT_FALSE(Polygon2d::compute_convex_hull({}, &polygon));
    EXPECT_FALSE(Polygon2d::compute_convex_hull({{1, 2}}, &polygon));
    EXPECT_FALSE(Polygon2d::compute_convex_hull({{3, 4}, {5, 6}}, &polygon));
    EXPECT_FALSE(
        Polygon2d::compute_convex_hull({{3, 4}, {3, 4}, {5, 6}, {5, 6}}, &polygon));

    EXPECT_TRUE(Polygon2d::compute_convex_hull({{0, 0}, {0, 4}, {3, 0}}, &polygon));
    EXPECT_TRUE(polygon.is_convex());
    EXPECT_NEAR(6.0, polygon.area(), 1e-5);

    EXPECT_TRUE(
        Polygon2d::compute_convex_hull({{0, 0}, {0, 4}, {3, 0}, {3, 4}}, &polygon));
    EXPECT_TRUE(polygon.is_convex());
    EXPECT_NEAR(12.0, polygon.area(), 1e-5);

    EXPECT_TRUE(Polygon2d::compute_convex_hull(
        {{0, 0}, {2, 2}, {1, 1}, {0, 4}, {3, 0}, {3, 4}}, &polygon));
    EXPECT_TRUE(polygon.is_convex());
    EXPECT_NEAR(12.0, polygon.area(), 1e-5);

    EXPECT_TRUE(Polygon2d::compute_convex_hull(
        {{0, 0}, {0, 4}, {0, 1}, {0, 3}, {0, 2}, {1, 0}, {3, 0}, {2, 0}},
        &polygon));
    EXPECT_TRUE(polygon.is_convex());
    EXPECT_NEAR(6.0, polygon.area(), 1e-5);

    for (int iter = 0; iter < 10000; ++iter) {
      const int kRange = 10;
      const int n = random_int(3, 10);
      std::vector<Vec2d> points;
      for (int i = 0; i < n; ++i) {
        points.emplace_back(random_int(0, kRange), random_int(0, kRange));
      }
      double area = 0;
      for (int x0 = 0; x0 < kRange; ++x0) {
        double min_y = 0.0;
        double max_y = 0.0;
        if (project_by_x_slow(points, x0 + 0.5, &min_y, &max_y)) {
          area += max_y - min_y;
        }
      }
      Polygon2d polygon;
      if (area < 1e-3) {
        EXPECT_FALSE(Polygon2d::compute_convex_hull(points, &polygon));
      } else {
        EXPECT_TRUE(Polygon2d::compute_convex_hull(points, &polygon));
        EXPECT_NEAR(area, polygon.area(), 1e-5);
      }
    }
  }
  TEST_END("convex hull");

  TEST_START("overlap");
  {
    const Polygon2d poly1(Box2d::create_aabox({0, 0}, {2, 2}));
    const Polygon2d poly2(Box2d::create_aabox({1, 1}, {3, 3}));
    const Polygon2d poly3(Box2d::create_aabox({2, 0}, {4, 2}));
    const Polygon2d poly4(Box2d({2, 2}, M_PI_4, sqrt(2.0), sqrt(2.0)));
    Polygon2d overlap_polygon;

    EXPECT_TRUE(poly1.compute_overlap(poly2, &overlap_polygon));
    EXPECT_NEAR(overlap_polygon.area(), 1.0, 1e-5);
    EXPECT_TRUE(poly2.compute_overlap(poly1, &overlap_polygon));
    EXPECT_NEAR(overlap_polygon.area(), 1.0, 1e-5);

    EXPECT_TRUE(poly2.compute_overlap(poly3, &overlap_polygon));
    EXPECT_NEAR(overlap_polygon.area(), 1.0, 1e-5);
    EXPECT_TRUE(poly3.compute_overlap(poly2, &overlap_polygon));
    EXPECT_NEAR(overlap_polygon.area(), 1.0, 1e-5);

    EXPECT_FALSE(poly1.compute_overlap(poly3, &overlap_polygon));
    EXPECT_FALSE(poly3.compute_overlap(poly1, &overlap_polygon));

    EXPECT_TRUE(poly1.compute_overlap(poly4, &overlap_polygon));
    EXPECT_NEAR(overlap_polygon.area(), 0.5, 1e-5);
    EXPECT_TRUE(poly4.compute_overlap(poly1, &overlap_polygon));
    EXPECT_NEAR(overlap_polygon.area(), 0.5, 1e-5);

    EXPECT_TRUE(poly2.compute_overlap(poly4, &overlap_polygon));
    EXPECT_NEAR(overlap_polygon.area(), 2.0, 1e-5);
    EXPECT_TRUE(poly4.compute_overlap(poly2, &overlap_polygon));
    EXPECT_NEAR(overlap_polygon.area(), 2.0, 1e-5);

    EXPECT_TRUE(poly3.compute_overlap(poly4, &overlap_polygon));
    EXPECT_NEAR(overlap_polygon.area(), 0.5, 1e-5);
    EXPECT_TRUE(poly4.compute_overlap(poly3, &overlap_polygon));
    EXPECT_NEAR(overlap_polygon.area(), 0.5, 1e-5);

    Vec2d first_intersect;
    Vec2d last_intersect;
    EXPECT_FALSE(poly1.get_overlap(LineSegment2d({-1, 0}, {-1, 2}),
                                  &first_intersect, &last_intersect));
    EXPECT_FALSE(poly1.get_overlap(LineSegment2d({-1, 1}, {-3, 1}),
                                  &first_intersect, &last_intersect));
    EXPECT_FALSE(poly1.get_overlap(LineSegment2d({1, 3}, {1, 5}), &first_intersect,
                                  &last_intersect));

    EXPECT_TRUE(poly1.get_overlap(LineSegment2d({1, -1}, {1, 3}), &first_intersect,
                                &last_intersect));
    EXPECT_NEAR(1.0, first_intersect.x(), 1e-5);
    EXPECT_NEAR(0.0, first_intersect.y(), 1e-5);
    EXPECT_NEAR(1.0, last_intersect.x(), 1e-5);
    EXPECT_NEAR(2.0, last_intersect.y(), 1e-5);

    EXPECT_TRUE(poly1.get_overlap(LineSegment2d({1, 1}, {1, 3}), &first_intersect,
                                &last_intersect));
    EXPECT_NEAR(1.0, first_intersect.x(), 1e-5);
    EXPECT_NEAR(1.0, first_intersect.y(), 1e-5);
    EXPECT_NEAR(1.0, last_intersect.x(), 1e-5);
    EXPECT_NEAR(2.0, last_intersect.y(), 1e-5);

    EXPECT_TRUE(poly1.get_overlap(LineSegment2d({1, -1}, {1, 1}), &first_intersect,
                                &last_intersect));
    EXPECT_NEAR(1.0, first_intersect.x(), 1e-5);
    EXPECT_NEAR(0.0, first_intersect.y(), 1e-5);
    EXPECT_NEAR(1.0, last_intersect.x(), 1e-5);
    EXPECT_NEAR(1.0, last_intersect.y(), 1e-5);

    EXPECT_TRUE(poly1.get_overlap(LineSegment2d({1, 3}, {3, 1}), &first_intersect,
                                &last_intersect));
    EXPECT_NEAR(2.0, first_intersect.x(), 1e-5);
    EXPECT_NEAR(2.0, first_intersect.y(), 1e-5);
    EXPECT_NEAR(2.0, last_intersect.x(), 1e-5);
    EXPECT_NEAR(2.0, last_intersect.y(), 1e-5);

    EXPECT_FALSE(poly1.get_overlap(LineSegment2d({4, 3}, {4, 3}), &first_intersect,
                                  &last_intersect));
    EXPECT_TRUE(poly1.get_overlap(LineSegment2d({1, 1}, {1, 1}), &first_intersect,
                                &last_intersect));
    EXPECT_NEAR(1.0, first_intersect.x(), 1e-5);
    EXPECT_NEAR(1.0, first_intersect.y(), 1e-5);
    EXPECT_NEAR(1.0, last_intersect.x(), 1e-5);
    EXPECT_NEAR(1.0, last_intersect.y(), 1e-5);

    const Polygon2d poly5(
        {{0, 0}, {4, 0}, {4, 2}, {3, 2}, {2, 1}, {1, 2}, {0, 2}});
    std::vector<LineSegment2d> overlap_line_segments =
        poly5.get_all_overlaps(LineSegment2d({-10, 1.5}, {10, 1.5}));
    EXPECT_EQ(2, overlap_line_segments.size());
    EXPECT_NEAR(0.0, overlap_line_segments[0].start().x(), 1e-5);
    EXPECT_NEAR(1.5, overlap_line_segments[0].start().y(), 1e-5);
    EXPECT_NEAR(1.5, overlap_line_segments[0].end().x(), 1e-5);
    EXPECT_NEAR(1.5, overlap_line_segments[0].end().y(), 1e-5);
    EXPECT_NEAR(2.5, overlap_line_segments[1].start().x(), 1e-5);
    EXPECT_NEAR(1.5, overlap_line_segments[1].start().y(), 1e-5);
    EXPECT_NEAR(4.0, overlap_line_segments[1].end().x(), 1e-5);
    EXPECT_NEAR(1.5, overlap_line_segments[1].end().y(), 1e-5);

    overlap_line_segments =
        poly5.get_all_overlaps(LineSegment2d({-10, 1}, {10, 1}));
    EXPECT_EQ(1, overlap_line_segments.size());
    EXPECT_NEAR(0.0, overlap_line_segments[0].start().x(), 1e-5);
    EXPECT_NEAR(1.0, overlap_line_segments[0].start().y(), 1e-5);
    EXPECT_NEAR(4.0, overlap_line_segments[0].end().x(), 1e-5);
    EXPECT_NEAR(1.0, overlap_line_segments[0].end().y(), 1e-5);

    overlap_line_segments =
        poly5.get_all_overlaps(LineSegment2d({-10, 0.5}, {10, 0.5}));
    EXPECT_EQ(1, overlap_line_segments.size());
    EXPECT_NEAR(0.0, overlap_line_segments[0].start().x(), 1e-5);
    EXPECT_NEAR(0.5, overlap_line_segments[0].start().y(), 1e-5);
    EXPECT_NEAR(4.0, overlap_line_segments[0].end().x(), 1e-5);
    EXPECT_NEAR(0.5, overlap_line_segments[0].end().y(), 1e-5);

    overlap_line_segments =
        poly5.get_all_overlaps(LineSegment2d({-10, -0.5}, {10, -0.5}));
    EXPECT_EQ(0, overlap_line_segments.size());
    overlap_line_segments =
        poly5.get_all_overlaps(LineSegment2d({-10, 2.5}, {10, 2.5}));
    EXPECT_EQ(0, overlap_line_segments.size());

    overlap_line_segments =
        poly5.get_all_overlaps(LineSegment2d({2, 0.5}, {2, 0.5}));
    EXPECT_EQ(1, overlap_line_segments.size());
    EXPECT_NEAR(2.0, overlap_line_segments[0].start().x(), 1e-5);
    EXPECT_NEAR(0.5, overlap_line_segments[0].start().y(), 1e-5);
    EXPECT_NEAR(2.0, overlap_line_segments[0].end().x(), 1e-5);
    EXPECT_NEAR(0.5, overlap_line_segments[0].end().y(), 1e-5);
    overlap_line_segments =
        poly5.get_all_overlaps(LineSegment2d({5, 0.5}, {5, 0.5}));
    EXPECT_EQ(0, overlap_line_segments.size());

    EXPECT_TRUE(Polygon2d({{0, 0}, {0, 4}, {4, 0}})
                    .compute_overlap(Polygon2d({{1, 1}, {1, 3}, {3, 1}}),
                                    &overlap_polygon));
    EXPECT_NEAR(overlap_polygon.area(), 2.0, 1e-5);

    EXPECT_TRUE(Polygon2d({{0, 0}, {0, 4}, {4, 0}})
                    .compute_overlap(Polygon2d({{1, 1}, {-1, 1}, {1, 3}}),
                                    &overlap_polygon));
    EXPECT_NEAR(overlap_polygon.area(), 1.5, 1e-5);
    EXPECT_TRUE(Polygon2d({{0, 0}, {0, 4}, {4, 0}})
                    .compute_overlap(Polygon2d({{2, 1}, {-1, 1}, {2, 4}}),
                                    &overlap_polygon));
    EXPECT_NEAR(overlap_polygon.area(), 3.0, 1e-5);
    EXPECT_TRUE(Polygon2d({{0, 0}, {0, 4}, {4, 0}})
                    .compute_overlap(Polygon2d({{3, 1}, {-1, 1}, {3, 5}}),
                                    &overlap_polygon));
    EXPECT_NEAR(overlap_polygon.area(), 3.5, 1e-5);
    EXPECT_TRUE(Polygon2d({{0, 0}, {0, 4}, {4, 0}})
                    .compute_overlap(Polygon2d({{4, 1}, {-1, 1}, {4, 6}}),
                                    &overlap_polygon));
    EXPECT_NEAR(overlap_polygon.area(), 3.5, 1e-5);

    for (int iter = 0; iter < 10000; ++iter) {
      const double x1 = random_double(-10, 10);
      const double y1 = random_double(-10, 10);
      const double x2 = random_double(-10, 10);
      const double y2 = random_double(-10, 10);
      const double heading1 = random_double(0, M_PI * 2.0);
      const double heading2 = random_double(0, M_PI * 2.0);
      const double l1 = random_double(1, 5);
      const double l2 = random_double(1, 5);
      const double w1 = random_double(1, 5);
      const double w2 = random_double(1, 5);
      Box2d box1({x1, y1}, heading1, l1, w1);
      Box2d box2({x2, y2}, heading2, l2, w2);
      Box2d shrinked_box2({x2, y2}, heading2, l2 - 0.2, w2 - 0.2);
      Box2d extended_box2({x2, y2}, heading2, l2 + 0.2, w2 + 0.2);
      if (!box1.has_overlap(extended_box2)) {
        EXPECT_FALSE(
            Polygon2d(box1).compute_overlap(Polygon2d(box2), &overlap_polygon));
      } else if (box1.has_overlap(shrinked_box2)) {
        EXPECT_TRUE(
            Polygon2d(box1).compute_overlap(Polygon2d(box2), &overlap_polygon));
      }
    }

    for (int iter = 0; iter < 10000; ++iter) {
      const int kRange = 10;
      const int n1 = random_int(3, 10);
      const int n2 = random_int(3, 10);
      std::vector<Vec2d> points1, points2;
      for (int i = 0; i < n1; ++i) {
        points1.emplace_back(random_int(0, kRange), random_int(0, kRange));
      }
      for (int i = 0; i < n2; ++i) {
        points2.emplace_back(random_int(0, kRange), random_int(0, kRange));
      }
      Polygon2d polygon1;
      Polygon2d polygon2;
      if (!Polygon2d::compute_convex_hull(points1, &polygon1) ||
          !Polygon2d::compute_convex_hull(points2, &polygon2)) {
        continue;
      }
      std::vector<double> key_points;
      for (int x0 = 0; x0 <= kRange; ++x0) {
        key_points.push_back(x0);
      }
      for (const auto &line_segment1 : polygon1.line_segments()) {
        for (const auto &line_segment2 : polygon2.line_segments()) {
          Vec2d pt;
          if (line_segment1.get_intersect(line_segment2, &pt)) {
            key_points.push_back(pt.x());
          }
        }
      }
      double area = 0;
      std::sort(key_points.begin(), key_points.end());
      for (size_t i = 0; i + 1 < key_points.size(); ++i) {
        const double width = key_points[i + 1] - key_points[i];
        if (width < 1e-6) {
          continue;
        }
        const double x = (key_points[i] + key_points[i + 1]) / 2.0;
        double min_y1 = 0.0;
        double max_y1 = 0.0;
        double min_y2 = 0.0;
        double max_y2 = 0.0;
        if (project_by_x_slow(points1, x, &min_y1, &max_y1) &&
            project_by_x_slow(points2, x, &min_y2, &max_y2)) {
          area +=
              std::max(0.0, std::min(max_y1, max_y2) - std::max(min_y1, min_y2)) *
              width;
        }
      }
      Polygon2d overlap_polygon;
      if (area < 1e-3) {
        EXPECT_FALSE(polygon1.compute_overlap(polygon2, &overlap_polygon));
      } else {
        EXPECT_TRUE(Polygon2d::compute_convex_hull(points1, &polygon1));
        EXPECT_TRUE(Polygon2d::compute_convex_hull(points2, &polygon2));
        EXPECT_TRUE(polygon1.compute_overlap(polygon2, &overlap_polygon));
        EXPECT_NEAR(area, overlap_polygon.area(), 1e-5);
      }
    }
  }
  TEST_END("overlap");

  TEST_START("bounding box");
  {
    Polygon2d poly1(Box2d::create_aabox({0, 0}, {2, 2}));
    Box2d box = poly1.bounding_box_with_heading(0.0);
    EXPECT_NEAR(1.0, box.center().x(), 1e-5);
    EXPECT_NEAR(1.0, box.center().y(), 1e-5);
    EXPECT_NEAR(4.0, box.area(), 1e-5);
    EXPECT_TRUE(Polygon2d(box).contains(poly1));
    AABox2d aabox = poly1.aabounding_box();
    EXPECT_NEAR(1.0, aabox.center().x(), 1e-5);
    EXPECT_NEAR(1.0, aabox.center().y(), 1e-5);
    EXPECT_NEAR(4.0, aabox.area(), 1e-5);
    EXPECT_NEAR(2.0, aabox.length(), 1e-5);
    EXPECT_NEAR(2.0, aabox.width(), 1e-5);

    box = poly1.bounding_box_with_heading(M_PI_4);
    EXPECT_NEAR(1.0, box.center().x(), 1e-5);
    EXPECT_NEAR(1.0, box.center().y(), 1e-5);
    EXPECT_NEAR(8.0, box.area(), 1e-5);
    EXPECT_TRUE(Polygon2d(box).contains(poly1));

    box = poly1.min_area_bounding_box();
    EXPECT_NEAR(1.0, box.center().x(), 1e-5);
    EXPECT_NEAR(1.0, box.center().y(), 1e-5);
    EXPECT_NEAR(4.0, box.area(), 1e-5);
    EXPECT_TRUE(Polygon2d(box).contains(poly1));

    Polygon2d poly2({{1, 0}, {0, 1}, {-1, 0}, {0, -1}});
    box = poly2.bounding_box_with_heading(0.0);
    EXPECT_NEAR(0.0, box.center().x(), 1e-5);
    EXPECT_NEAR(0.0, box.center().y(), 1e-5);
    EXPECT_NEAR(4.0, box.area(), 1e-5);
    EXPECT_TRUE(Polygon2d(box).contains(poly2));
    aabox = poly2.aabounding_box();
    EXPECT_NEAR(0.0, aabox.center().x(), 1e-5);
    EXPECT_NEAR(0.0, aabox.center().y(), 1e-5);
    EXPECT_NEAR(4.0, aabox.area(), 1e-5);
    EXPECT_NEAR(2.0, aabox.length(), 1e-5);
    EXPECT_NEAR(2.0, aabox.width(), 1e-5);

    box = poly2.bounding_box_with_heading(M_PI_4);
    EXPECT_NEAR(0.0, box.center().x(), 1e-5);
    EXPECT_NEAR(0.0, box.center().y(), 1e-5);
    EXPECT_NEAR(2.0, box.area(), 1e-5);
    EXPECT_TRUE(Polygon2d(box).contains(poly2));

    box = poly2.min_area_bounding_box();
    EXPECT_NEAR(0.0, box.center().x(), 1e-5);
    EXPECT_NEAR(0.0, box.center().y(), 1e-5);
    EXPECT_NEAR(2.0, box.area(), 1e-5);
    EXPECT_TRUE(Polygon2d(box).contains(poly2));

    for (int iter = 0; iter < 1000; ++iter) {
      const int num_sample_points = random_int(3, 10);
      std::vector<Vec2d> points;
      for (int i = 0; i < num_sample_points; ++i) {
        const double x = random_int(-10, 10);
        const double y = random_int(-10, 10);
        points.emplace_back(x, y);
      }
      Polygon2d polygon;
      if (!Polygon2d::compute_convex_hull(points, &polygon)) {
        continue;
      }
      double min_area = std::numeric_limits<double>::infinity();
      for (int iter2 = 0; iter2 < 10; ++iter2) {
        const double heading = random_double(0, M_PI * 2.0);
        box = polygon.bounding_box_with_heading(heading);
        EXPECT_TRUE(Polygon2d(box).contains(polygon));
        min_area = std::min(min_area, box.area());
      }
      box = polygon.min_area_bounding_box();
      EXPECT_TRUE(Polygon2d(box).contains(polygon));
      EXPECT_LE(box.area(), min_area + 1e-5);
    }
  }
  TEST_END("bounding box");

  TEST_START("expand");
  {
    {
      const Polygon2d poly(Box2d::create_aabox({0, 0}, {2, 2}));
      const Polygon2d exp_poly = poly.expand_by_distance(1.0);
      EXPECT_TRUE(exp_poly.is_convex());
      const Box2d box = exp_poly.bounding_box_with_heading(0.0);
      EXPECT_NEAR(box.center().x(), 1.0, 1e-6);
      EXPECT_NEAR(box.center().y(), 1.0, 1e-6);
      EXPECT_NEAR(box.width(), 4.0, 1e-6);
      EXPECT_NEAR(box.length(), 4.0, 1e-6);
      EXPECT_NEAR(exp_poly.area(), 12 + M_PI, 0.1);
    }
    {
      const std::vector<Vec2d> points{{0, 0}, {2, 0}, {2, 2}, {0, 2}, {1, 1}};
      const Polygon2d poly(points);
      const Polygon2d exp_poly = poly.expand_by_distance(1.0);
      EXPECT_TRUE(exp_poly.is_convex());
      const Box2d box = exp_poly.bounding_box_with_heading(0.0);
      EXPECT_NEAR(box.center().x(), 1.0, 1e-6);
      EXPECT_NEAR(box.center().y(), 1.0, 1e-6);
      EXPECT_NEAR(box.width(), 4.0, 1e-6);
      EXPECT_NEAR(box.length(), 4.0, 1e-6);
      EXPECT_NEAR(exp_poly.area(), 12 + M_PI, 0.1);
    }
  }
  TEST_END("expand");
}
