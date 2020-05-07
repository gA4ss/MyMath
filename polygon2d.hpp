#ifndef MYMATH_POLYGON2D_HPP
#define MYMATH_POLYGON2D_HPP

#include "mymath_config.h"

#include <cmath>
#include <cassert>
#include <vector>
#include <limits>
#include <utility>
#include <algorithm>

#ifdef MYMATH_DBG
#include <string>
#include <sstream>
#endif

#include "box2d.hpp"
#include "line_segment2d.hpp"
#include "vec2d.hpp"
#include "math_utils.hpp"

namespace mypilot {
namespace mymath {

// 2D-多边形
class Polygon2d {
public:
  Polygon2d() = default;

  explicit Polygon2d(const Box2d& box) {
    box.get_all_corners(&_points);
    build_from_points();
  }

  explicit Polygon2d(std::vector<Vec2d> points) : _points(std::move(points)) {
    build_from_points();
  }

  const std::vector<Vec2d>& points() const { return _points; }
  const std::vector<LineSegment2d>& line_segments() const { return _line_segments; }
  int num_points() const { return _num_points; }
  bool is_convex() const { return _is_convex; } // 多边形是否是凸的
  double area() const { return _area; }

  // 计算点到多边形最短距离,如果点在多边形内返回0。
  double distance_to(const Vec2d& point) const {
    assert(_points.size() >= 3);
    if (is_point_in(point)) {
      return 0.0;
    }
    double distance = std::numeric_limits<double>::infinity();
    for (int i = 0; i < _num_points; ++i) {
      distance = std::min(distance, _line_segments[i].distance_to(point));
    }
    return distance;
  }

  // 计算点到多边形最短距离的平方,如果点在多边形内返回0。
  double distance_square_to(const Vec2d& point) const {
    assert(_points.size() >= 3);
    if (is_point_in(point)) {
      return 0.0;
    }
    double distance_sqr = std::numeric_limits<double>::infinity();
    for (int i = 0; i < _num_points; ++i) {
      distance_sqr =
        std::min(distance_sqr, _line_segments[i].distance_square_to(point));
    }
    return distance_sqr;
  }

  // 计算线段到多边形最短距离,如果线段与多边形存在交集返回0。
  double distance_to(const LineSegment2d& line_segment) const {
    if (line_segment.length() <= math_epsilon) {
      return distance_to(line_segment.start());
    }
    assert(_points.size() >= 3);
    if (is_point_in(line_segment.center())) {
      return 0.0;
    }
    if (std::any_of(_line_segments.begin(), _line_segments.end(),
                    [&](const LineSegment2d& poly_seg) {
                      return poly_seg.has_intersect(line_segment);
                    })) {
      return 0.0;
    }

    double distance = std::min(distance_to(line_segment.start()),
                               distance_to(line_segment.end()));
    for (int i = 0; i < _num_points; ++i) {
      distance = std::min(distance, line_segment.distance_to(_points[i]));
    }
    return distance;
  }

  // 计算一个方形盒子到多边形最短距离,如果方形盒子与多边形存在交集返回0。
  double distance_to(const Box2d& box) const {
    assert(_points.size() >= 3);
    return distance_to(Polygon2d(box));
  }

  // 计算一个多边形到多边形最短距离,如果多边形与多边形存在交集返回0。
  double distance_to(const Polygon2d& polygon) const {
    assert(_points.size() >= 3);
    assert(polygon.num_points() >= 3);

    if (is_point_in(polygon.points()[0])) {
      return 0.0;
    }
    if (polygon.is_point_in(_points[0])) {
      return 0.0;
    }
    double distance = std::numeric_limits<double>::infinity();
    for (int i = 0; i < _num_points; ++i) {
      distance = std::min(distance, polygon.distance_to(_line_segments[i]));
    }
    return distance;
  }

  // 计算一个点到多边形的最短距离
  double distance_to_boundary(const Vec2d& point) const {
    double distance = std::numeric_limits<double>::infinity();
    for (int i = 0; i < _num_points; ++i) {
      distance = std::min(distance, _line_segments[i].distance_to(point));
    }
    return distance;
  }

  // 给定点是否在多边形的边上
  bool is_point_on_boundary(const Vec2d& point) const {
    assert(_points.size() >= 3);
    return std::any_of(
        _line_segments.begin(), _line_segments.end(),
        [&](const LineSegment2d& poly_seg) { return poly_seg.is_point_in(point); });
  }

  // 计算点是否在多边形内
  bool is_point_in(const Vec2d& point) const {
    assert(_points.size() >= 3);
    if (is_point_on_boundary(point)) {
      return true;
    }
    int j = _num_points - 1;
    int c = 0;
    for (int i = 0; i < _num_points; ++i) {
      if ((_points[i].y() > point.y()) != (_points[j].y() > point.y())) {
        const double side = cross_prod(point, _points[i], _points[j]);
        if (_points[i].y() < _points[j].y() ? side > 0.0 : side < 0.0) {
          ++c;
        }
      }
      j = i;
    }
    return c & 1;
  }

  // 检查多边形是否包含一条线段
  bool contains(const LineSegment2d& line_segment) const {
    if (line_segment.length() <= math_epsilon) {
      return is_point_in(line_segment.start());
    }
    assert(_points.size() >= 3);
    if (!is_point_in(line_segment.start())) {
      return false;
    }
    if (!is_point_in(line_segment.end())) {
      return false;
    }
    if (!_is_convex) {
      std::vector<LineSegment2d> overlaps = get_all_overlaps(line_segment);
      double total_length = 0;
      for (const auto &overlap_seg : overlaps) {
        total_length += overlap_seg.length();
      }
      return total_length >= line_segment.length() - math_epsilon;
    }
    return true;
  }

  // 是否包含目标多边形
  bool contains(const Polygon2d& polygon) const {
    assert(_points.size() >= 3);
    if (_area < polygon.area() - math_epsilon) {
      return false;
    }
    if (!is_point_in(polygon.points()[0])) {
      return false;
    }
    const auto& line_segments = polygon.line_segments();
    return std::all_of(line_segments.begin(), line_segments.end(),
                      [&](const LineSegment2d &line_segment) {
                        return contains(line_segment);
                      });
  }

  // 是否与线段包含重回
  bool has_overlap(const LineSegment2d& line_segment) const {
    assert(_points.size() >= 3);
    if ((line_segment.start().x() < _min_x && line_segment.end().x() < _min_x) ||
        (line_segment.start().x() > _max_x && line_segment.end().x() > _max_x) ||
        (line_segment.start().y() < _min_y && line_segment.end().y() < _min_y) ||
        (line_segment.start().y() > _max_y && line_segment.end().y() > _max_y)) {
      return false;
    }
    Vec2d first;
    Vec2d last;
    return get_overlap(line_segment, &first, &last);
  }

  // 获取与线段的重合，返回重叠部分的端点。
  bool get_overlap(const LineSegment2d& line_segment, Vec2d* const first,
                   Vec2d* const last) const {
    assert(_points.size() >= 3);
    assert(first);
    assert(last);

    if (line_segment.length() <= math_epsilon) {
      if (!is_point_in(line_segment.start())) {
        return false;
      }
      *first = line_segment.start();
      *last = line_segment.start();
      return true;
    }

    double min_proj = line_segment.length();
    double max_proj = 0;
    if (is_point_in(line_segment.start())) {
      *first = line_segment.start();
      min_proj = 0.0;
    }
    if (is_point_in(line_segment.end())) {
      *last = line_segment.end();
      max_proj = line_segment.length();
    }
    for (const auto& poly_seg : _line_segments) {
      Vec2d pt;
      if (poly_seg.get_intersect(line_segment, &pt)) {
        const double proj = line_segment.project_onto_unit(pt);
        if (proj < min_proj) {
          min_proj = proj;
          *first = pt;
        }
        if (proj > max_proj) {
          max_proj = proj;
          *last = pt;
        }
      }
    }
    return min_proj <= max_proj + math_epsilon;
  }

  // 获取线段与多边形的重合,返回线段表示。
  std::vector<LineSegment2d> get_all_overlaps(const LineSegment2d& line_segment) const {
    assert(_points.size() >= 3);

    if (line_segment.length() <= math_epsilon) {
      std::vector<LineSegment2d> overlaps;
      if (is_point_in(line_segment.start())) {
        overlaps.push_back(line_segment);
      }
      return overlaps;
    }
    std::vector<double> projections;
    if (is_point_in(line_segment.start())) {
      projections.push_back(0.0);
    }
    if (is_point_in(line_segment.end())) {
      projections.push_back(line_segment.length());
    }
    for (const auto& poly_seg : _line_segments) {
      Vec2d pt;
      if (poly_seg.get_intersect(line_segment, &pt)) {
        projections.push_back(line_segment.project_onto_unit(pt));
      }
    }
    std::sort(projections.begin(), projections.end());
    std::vector<std::pair<double, double>> overlaps;
    for (size_t i = 0; i + 1 < projections.size(); ++i) {
      const double start_proj = projections[i];
      const double end_proj = projections[i + 1];
      if (end_proj - start_proj <= math_epsilon) {
        continue;
      }
      const Vec2d reference_point =
          line_segment.start() +
          (start_proj + end_proj) / 2.0 * line_segment.unit_direction();
      if (!is_point_in(reference_point)) {
        continue;
      }
      if (overlaps.empty() ||
          start_proj > overlaps.back().second + math_epsilon) {
        overlaps.emplace_back(start_proj, end_proj);
      } else {
        overlaps.back().second = end_proj;
      }
    }
    std::vector<LineSegment2d> overlap_line_segments;
    for (const auto &overlap : overlaps) {
      overlap_line_segments.emplace_back(
          line_segment.start() + overlap.first * line_segment.unit_direction(),
          line_segment.start() + overlap.second * line_segment.unit_direction());
    }
    return overlap_line_segments;
  }

  // 获取所有的顶点坐标
  void get_all_vertices(std::vector<Vec2d>* const vertices) const { if (vertices) *vertices = _points; }
  std::vector<Vec2d> get_all_vertices() const { return _points; }

  // 判断两个多边形是否有重叠
  bool has_overlap(const Polygon2d& polygon) const {
    assert(_points.size() >= 3);
    if (polygon.max_x() < min_x() || polygon.min_x() > max_x() ||
        polygon.max_y() < min_y() || polygon.min_y() > max_y()) {
      return false;
    }
    return distance_to(polygon) <= math_epsilon;
  }

  // 计算两个多边形的重叠并返回重叠部分(仅当两个多边形是凸多边形时才计算)
  bool compute_overlap(const Polygon2d& other_polygon,
                       Polygon2d* const overlap_polygon) const {
    assert(_points.size() >= 3);
    assert(overlap_polygon);
    assert(_is_convex && other_polygon.is_convex());
    std::vector<Vec2d> points = other_polygon.points();
    for (int i = 0; i < _num_points; ++i) {
      if (!clip_convex_hull(_line_segments[i], &points)) {
        return false;
      }
    }
    return compute_convex_hull(points, overlap_polygon);
  }

  // 按照当前多边形返回一个AABoundingBox
  AABox2d aabounding_box() const {
    return AABox2d({_min_x, _min_y}, {_max_x, _max_y});
  }

  // 根据偏航角'heading'返回一个长方形结构
  Box2d bounding_box_with_heading(const double heading) const {
    assert(_points.size() >= 3);
    const Vec2d direction_vec = Vec2d::create_unit_vec2d(heading);
    Vec2d px1;
    Vec2d px2;
    Vec2d py1;
    Vec2d py2;
    extreme_points(heading, &px1, &px2);
    extreme_points(heading - M_PI_2, &py1, &py2);
    const double x1 = px1.inner_prod(direction_vec);
    const double x2 = px2.inner_prod(direction_vec);
    const double y1 = py1.cross_prod(direction_vec);
    const double y2 = py2.cross_prod(direction_vec);
    return Box2d(
      (x1 + x2) / 2.0 * direction_vec +
        (y1 + y2) / 2.0 * Vec2d(direction_vec.y(), -direction_vec.x()),
      heading, x2 - x1, y2 - y1);
  }

  // 根据多边形返回最小面积的长方形
  Box2d min_area_bounding_box() const {
    assert(_points.size() >= 3);
    if (!_is_convex) {
      Polygon2d convex_polygon;
      compute_convex_hull(_points, &convex_polygon);
      assert(convex_polygon.is_convex());
      return convex_polygon.min_area_bounding_box();
    }
    double min_area = std::numeric_limits<double>::infinity();
    double min_area_at_heading = 0.0;
    int left_most = 0;
    int right_most = 0;
    int top_most = 0;
    // 遍历多边形所有点
    for (int i = 0; i < _num_points; ++i) {
      const auto& line_segment = _line_segments[i];
      double proj = 0.0;
      double min_proj = line_segment.project_onto_unit(_points[left_most]);
      while ((proj = line_segment.project_onto_unit(_points[prev(left_most)])) <
            min_proj) {
        min_proj = proj;
        left_most = prev(left_most);
      }
      while ((proj = line_segment.project_onto_unit(_points[next(left_most)])) <
            min_proj) {
        min_proj = proj;
        left_most = next(left_most);
      }
      double max_proj = line_segment.project_onto_unit(_points[right_most]);
      while ((proj = line_segment.project_onto_unit(_points[prev(right_most)])) >
            max_proj) {
        max_proj = proj;
        right_most = prev(right_most);
      }
      while ((proj = line_segment.project_onto_unit(_points[next(right_most)])) >
            max_proj) {
        max_proj = proj;
        right_most = next(right_most);
      }
      double prod = 0.0;
      double max_prod = line_segment.project_onto_unit(_points[top_most]);
      while ((prod = line_segment.project_onto_unit(_points[prev(top_most)])) >
            max_prod) {
        max_prod = prod;
        top_most = prev(top_most);
      }
      while ((prod = line_segment.project_onto_unit(_points[next(top_most)])) >
            max_prod) {
        max_prod = prod;
        top_most = next(top_most);
      }
      const double area = max_prod * (max_proj - min_proj);
      if (area < min_area) {
        min_area = area;
        min_area_at_heading = line_segment.heading();
      }
    }
    return bounding_box_with_heading(min_area_at_heading);
  }

  // 获取沿航向的极点
  void extreme_points(const double heading, 
                      Vec2d* const first, Vec2d* const last) const {
    assert(_points.size() >= 3);
    assert(first);
    assert(last);

    // 创建方向向量
    const Vec2d direction_vec = Vec2d::create_unit_vec2d(heading);
    double min_proj = std::numeric_limits<double>::infinity();
    double max_proj = -std::numeric_limits<double>::infinity();

    // 遍历所有点
    for (const auto& pt : _points) {
      const double proj = pt.inner_prod(direction_vec);
      if (proj < min_proj) {
        min_proj = proj;
        *first = pt;
      }
      if (proj > max_proj) {
        max_proj = proj;
        *last = pt;
      }
    }
  }

  // 将此多边形扩展一定距离'distance'
  Polygon2d expand_by_distance(const double distance) const {
    if (!_is_convex) {
      Polygon2d convex_polygon;
      compute_convex_hull(_points, &convex_polygon);
      assert(convex_polygon.is_convex());
      return convex_polygon.expand_by_distance(distance);
    }
    const double min_angle = 0.1;
    std::vector<Vec2d> points;
    for (int i = 0; i < _num_points; ++i) {
      const double start_angle = _line_segments[prev(i)].heading() - M_PI_2;
      const double end_angle = _line_segments[i].heading() - M_PI_2;
      const double diff = wrap_angle(end_angle - start_angle);
      if (diff <= math_epsilon) {
        points.push_back(_points[i] +
          Vec2d::create_unit_vec2d(start_angle) * distance);
      } else {
        const int count = static_cast<int>(diff / min_angle) + 1;
        for (int k = 0; k <= count; ++k) {
          const double angle = start_angle +
            diff * static_cast<double>(k) / static_cast<double>(count);
          points.push_back(_points[i] + Vec2d::create_unit_vec2d(angle) * distance);
        }
      }
    }
    Polygon2d new_polygon;
    assert(compute_convex_hull(points, &new_polygon));
    return new_polygon;
  }

#if MYMATH_DBG
  std::string str() const {
    std::ostringstream ss;
    std::ostringstream pints_str;
    for (int i = 0; i < _points.size(); i++) {
      if (i < _points.size()-1)
        pints_str << "( " << _points[i].x() << ", " << _points[i].y() << " )" << ",";
      else
        pints_str << "( " << _points[i].x() << ", " << _points[i].y() << " )";
    }

    ss << "polygon2d (  num_points = " << _num_points
       << "  points = ( " << pints_str.str() << " )  "
       << (_is_convex ? "convex" : "non-convex")
       << "  area = " << _area << " )";
    
    return ss.str();
  }
#endif

  // 计算一组点'points'的凸包'polygon'
  static bool compute_convex_hull(const std::vector<Vec2d>& points,
                                  Polygon2d* const polygon) {
    assert(polygon);
    const int n = points.size();
    if (n < 3) {
      return false;
    }
    // 初始化索引
    std::vector<int> sorted_indices(n);
    for (int i = 0; i < n; ++i) {
      sorted_indices[i] = i;
    }
    // 对所有点进行排序
    std::sort(sorted_indices.begin(), sorted_indices.end(),
              [&](const int idx1, const int idx2) {
                const Vec2d& pt1 = points[idx1];
                const Vec2d& pt2 = points[idx2];
                const double dx = pt1.x() - pt2.x();
                if (std::abs(dx) > math_epsilon) {
                  return dx < 0.0;
                }
                return pt1.y() < pt2.y();
              });
    int count = 0;
    std::vector<int> results;
    results.reserve(n);
    int last_count = 1;
    for (int i = 0; i < n + n; ++i) {
      if (i == n) {
        last_count = count;
      }
      const int idx = sorted_indices[(i < n) ? i : (n + n - 1 - i)];
      const Vec2d& pt = points[idx];
      while (count > last_count &&
        cross_prod(points[results[count - 2]], points[results[count - 1]],
        pt) <= math_epsilon) {
        results.pop_back();
        --count;
      }
      results.push_back(idx);
      ++count;
    }
    --count;
    if (count < 3) {
      return false;
    }
    std::vector<Vec2d> result_points;
    result_points.reserve(count);
    for (int i = 0; i < count; ++i) {
      result_points.push_back(points[results[i]]);
    }
    *polygon = Polygon2d(result_points);
    return true;
  }

  static bool clip_convex_hull(const LineSegment2d& line_segment,
                               std::vector<Vec2d>* const points) {
    if (line_segment.length() <= math_epsilon) {
      return true;
    }
    assert(points);
    const int n = points->size();
    if (n < 3) {
      return false;
    }
    std::vector<double> prod(n);
    std::vector<int> side(n);
    for (int i = 0; i < n; ++i) {
      prod[i] = cross_prod(line_segment.start(), line_segment.end(), (*points)[i]);
      if (std::abs(prod[i]) <= math_epsilon) {
        side[i] = 0;
      } else {
        side[i] = ((prod[i] < 0) ? -1 : 1);
      }
    }

    std::vector<Vec2d> new_points;
    for (int i = 0; i < n; ++i) {
      if (side[i] >= 0) {
        new_points.push_back((*points)[i]);
      }
      const int j = ((i == n - 1) ? 0 : (i + 1));
      if (side[i] * side[j] < 0) {
        const double ratio = prod[j] / (prod[j] - prod[i]);
        new_points.emplace_back(
          (*points)[i].x() * ratio + (*points)[j].x() * (1.0 - ratio),
          (*points)[i].y() * ratio + (*points)[j].y() * (1.0 - ratio));
      }
    }

    points->swap(new_points);
    return points->size() >= 3;
  }

  double min_x() const { return _min_x; }
  double max_x() const { return _max_x; }
  double min_y() const { return _min_y; }
  double max_y() const { return _max_y; }

protected:
  void build_from_points() {
    _num_points = _points.size();
    assert(_num_points >= 3);

    // 确认点的顺序为'ccw'顺序
    _area = 0.0;
    for (int i = 1; i < _num_points; ++i) {
      _area += cross_prod(_points[0], _points[i - 1], _points[i]);
    }
    if (_area < 0) {
      _area = -_area;
      std::reverse(_points.begin(), _points.end());
    }
    _area /= 2.0;
    assert(_area > math_epsilon);

    // 构造线段
    _line_segments.reserve(_num_points);
    for (int i = 0; i < _num_points; ++i) {
      _line_segments.emplace_back(_points[i], _points[next(i)]);
    }

    // 检查凸度
    _is_convex = true;
    for (int i = 0; i < _num_points; ++i) {
      if (cross_prod(_points[prev(i)], 
                     _points[i], 
                     _points[next(i)]) <= -math_epsilon) {
        _is_convex = false;
        break;
      }
    }

    // 计算 aabox
    _min_x = _points[0].x();
    _max_x = _points[0].x();
    _min_y = _points[0].y();
    _max_y = _points[0].y();
    for (const auto &point : _points) {
      _min_x = std::min(_min_x, point.x());
      _max_x = std::max(_max_x, point.x());
      _min_y = std::min(_min_y, point.y());
      _max_y = std::max(_max_y, point.y());
    }
  }

  int next(int at) const {
    return at >= _num_points - 1 ? 0 : at + 1;
  }

  int prev(int at) const {
    return at == 0 ? _num_points - 1 : at - 1;
  }

  std::vector<Vec2d> _points;
  int _num_points = 0;
  std::vector<LineSegment2d> _line_segments;
  bool _is_convex = false;
  double _area = 0.0;
  double _min_x = 0.0;
  double _max_x = 0.0;
  double _min_y = 0.0;
  double _max_y = 0.0;
};

}}

#endif
