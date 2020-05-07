#include "aaboxkdtree2d.hpp"
#include "ltest.hpp"

#include <string>
#include <set>

#include "line_segment2d.hpp"
#include "math_utils.hpp"

using namespace mypilot::mymath;

class Object {
public:
  Object(const double x1, const double y1, const double x2, const double y2,
         const int id)
      : aabox_({x1, y1}, {x2, y2}),
        line_segment_({x1, y1}, {x2, y2}),
        id_(id) {}
  const AABox2d &aabox() const { return aabox_; }
  double distance_to(const Vec2d &point) const {
    return line_segment_.distance_to(point);
  }
  double distance_square_to(const Vec2d &point) const {
    return line_segment_.distance_square_to(point);
  }
  int id() const { return id_; }

private:
  AABox2d aabox_;
  LineSegment2d line_segment_;
  int id_ = 0;
};

int main(int argc, char* argv[]) {
  TEST_START("over all tests");
  {
    const int kNumBoxes[4] = {1, 10, 50, 100};
    const int kNumQueries = 1000;
    const double kSize = 100;
    const int kNumTrees = 4;
    AABoxKDTreeParams kdtree_params[kNumTrees];
    kdtree_params[1].max_depth = 2;
    kdtree_params[2].max_leaf_dimension = kSize / 4.0;
    kdtree_params[3].max_leaf_size = 20;

    for (int num_boxes : kNumBoxes) {
      std::vector<Object> objects;
      for (int i = 0; i < num_boxes; ++i) {
        const double cx = random_double(-kSize, kSize);
        const double cy = random_double(-kSize, kSize);
        const double dx = random_double(-kSize / 10.0, kSize / 10.0);
        const double dy = random_double(-kSize / 10.0, kSize / 10.0);
        objects.emplace_back(cx - dx, cy - dy, cx + dx, cy + dy, i);
      }
      std::unique_ptr<AABoxKDTree2d<Object>> kdtrees[kNumTrees];
      for (int i = 0; i < kNumTrees; ++i) {
        kdtrees[i].reset(new AABoxKDTree2d<Object>(objects, kdtree_params[i]));
      }
      for (int i = 0; i < kNumQueries; ++i) {
        const Vec2d point(random_double(-kSize * 1.5, kSize * 1.5),
                          random_double(-kSize * 1.5, kSize * 1.5));
        double expected_distance = std::numeric_limits<double>::infinity();
        for (const auto &object : objects) {
          expected_distance =
              std::min(expected_distance, object.distance_to(point));
        }
        for (int k = 0; k < kNumTrees; ++k) {
          const Object *nearest_object = kdtrees[k]->get_nearest_object(point);
          const double actual_distance = nearest_object->distance_to(point);
          EXPECT_NEAR(actual_distance, expected_distance, 1e-3);
        }
      }
      for (int i = 0; i < kNumQueries; ++i) {
        const Vec2d point(random_double(-kSize * 1.5, kSize * 1.5),
                          random_double(-kSize * 1.5, kSize * 1.5));
        const double distance = random_double(0, kSize * 2.0);
        for (int k = 0; k < kNumTrees; ++k) {
          std::vector<const Object *> result_objects =
              kdtrees[k]->get_objects(point, distance);
          std::set<int> result_ids;
          for (const Object *object : result_objects) {
            result_ids.insert(object->id());
          }
          EXPECT_EQ(result_objects.size(), result_ids.size());
          for (const auto &object : objects) {
            const double d = object.distance_to(point);
            if (std::abs(d - distance) <= 1e-3) {
              continue;
            }
            if (d < distance) {
              EXPECT_TRUE(result_ids.count(object.id()));
            } else {
              EXPECT_FALSE(result_ids.count(object.id()));
            }
          }
        }
      }
    }
  }
  TEST_END("over all tests");
}
