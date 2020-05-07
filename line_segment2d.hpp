#ifndef MYMATH_LINE_SEGEMENT2D_HPP
#define MYMATH_LINE_SEGEMENT2D_HPP

#include "mymath_config.h"

#include <cmath>
#include <cassert>
#include <utility>
#include <algorithm>

#ifdef MYMATH_DBG
#include <string>
#include <sstream>
#endif

#include "vec2d.hpp"
#include "math_utils.hpp"

namespace mypilot {
namespace mymath {

class LineSegment2d {
public:
  LineSegment2d() { _unit_direction = Vec2d(1, 0); }
  LineSegment2d(const Vec2d& start, const Vec2d& end) : 
    _start(start), _end(end) {
    const double dx = _end.x() - _start.x();
    const double dy = _end.y() - _start.y();
    _length = hypot(dx, dy);
    // 单位向量
    _unit_direction = (_length <= math_epsilon ? Vec2d(0, 0)
      : Vec2d(dx / _length, dy / _length));
    // 获取线段的角度
    _heading = _unit_direction.angle();
  }

  const Vec2d& start() const { return _start; }
  const Vec2d& end() const { return _end; }
  const Vec2d& unit_direction() const { return _unit_direction; }
  Vec2d center() const { return (_start + _end) / 2.0; }
  double heading() const { return _heading; }
  double cos_heading() const { return _unit_direction.x(); }
  double sin_heading() const { return _unit_direction.y(); }
  double length() const { return _length; }
  double length_sqr() const { return _length * _length; }

  // 计算当前线段与点的距离
  double distance_to(const Vec2d& point) const {
    // 如果当前线段长度非常小则当作点对待
    if (_length <= math_epsilon) {
      return point.distance_to(_start);
    }

    // 计算投影距离
    const double x0 = point.x() - _start.x();
    const double y0 = point.y() - _start.y();
    const double proj = x0 * _unit_direction.x() + y0 * _unit_direction.y();
    if (proj <= 0.0) {
      return hypot(x0, y0);
    }
    if (proj >= _length) {
      return point.distance_to(_end);
    }
    return std::abs(x0 * _unit_direction.y() - y0 * _unit_direction.x());
  }

  // 计算从线段上的点到2-D点的最短距离，并获得线段上最近的点。
  double distance_to(const Vec2d& point, Vec2d* const nearest_pt) const {
    assert(nearest_pt);
    if (_length <= math_epsilon) {
      *nearest_pt = _start;
      return point.distance_to(_start);
    }
    const double x0 = point.x() - _start.x();
    const double y0 = point.y() - _start.y();
    const double proj = x0 * _unit_direction.x() + y0 * _unit_direction.y();
    if (proj < 0.0) {
      *nearest_pt = _start;
      return hypot(x0, y0);
    }
    if (proj > _length) {
      *nearest_pt = _end;
      return point.distance_to(_end);
    }
    *nearest_pt = _start + _unit_direction * proj;
    return std::abs(x0 * _unit_direction.y() - y0 * _unit_direction.x());
  }

  // 计算从线段上的点到2-D点的最短距离的平方
  double distance_square_to(const Vec2d& point) const {
    if (_length <= math_epsilon) {
      return point.distance_square_to(_start);
    }
    const double x0 = point.x() - _start.x();
    const double y0 = point.y() - _start.y();
    const double proj = x0 * _unit_direction.x() + y0 * _unit_direction.y();
    if (proj <= 0.0) {
      return square(x0) + square(y0);
    }
    if (proj >= _length) {
      return point.distance_square_to(_end);
    }
    return square(x0 * _unit_direction.y() - y0 * _unit_direction.x());
  }

  // 计算从线段上的点到2-D点的最短距离的平方，并获得线段上最近的点。
  double distance_square_to(const Vec2d& point, Vec2d* const nearest_pt) const {
    assert(nearest_pt);
    if (_length <= math_epsilon) {
      *nearest_pt = _start;
      return point.distance_square_to(_start);
    }
    const double x0 = point.x() - _start.x();
    const double y0 = point.y() - _start.y();
    const double proj = x0 * _unit_direction.x() + y0 * _unit_direction.y();
    if (proj <= 0.0) {
      *nearest_pt = _start;
      return square(x0) + square(y0);
    }
    if (proj >= _length) {
      *nearest_pt = _end;
      return point.distance_square_to(_end);
    }
    *nearest_pt = _start + _unit_direction * proj;
    return square(x0 * _unit_direction.y() - y0 * _unit_direction.x());
  }

  // 确定点是否在线段上
  bool is_point_in(const Vec2d& point) const {
    if (_length <= math_epsilon) {
      return std::abs(point.x() - _start.x()) <= math_epsilon &&
        std::abs(point.y() - _start.y()) <= math_epsilon;
    }
    const double prod = cross_prod(point, _start, _end);
    if (std::abs(prod) > math_epsilon) {
      return false;
    }
    return is_with_in(point.x(), _start.x(), _end.x()) &&
      is_with_in(point.y(), _start.y(), _end.y());
  }

  // 检查两条线段是否有交点
  bool has_intersect(const LineSegment2d& other_segment) const {
    Vec2d point;
    return get_intersect(other_segment, &point);
  }

  // 过去与另外一条线段的交点
  bool get_intersect(const LineSegment2d& other_segment, 
                     Vec2d* const point) const {
    assert(point);
    if (is_point_in(other_segment.start())) {
      *point = other_segment.start();
      return true;
    }
    if (is_point_in(other_segment.end())) {
      *point = other_segment.end();
      return true;
    }
    if (other_segment.is_point_in(_start)) {
      *point = _start;
      return true;
    }
    if (other_segment.is_point_in(_end)) {
      *point = _end;
      return true;
    }
    if (_length <= math_epsilon || other_segment.length() <= math_epsilon) {
      return false;
    }
    const double cc1 = cross_prod(_start, _end, other_segment.start());
    const double cc2 = cross_prod(_start, _end, other_segment.end());
    if (cc1 * cc2 >= -math_epsilon) {
      return false;
    }
    const double cc3 =
      cross_prod(other_segment.start(), other_segment.end(), _start);
    const double cc4 =
      cross_prod(other_segment.start(), other_segment.end(), _end);
    if (cc3 * cc4 >= -math_epsilon) {
      return false;
    }
    const double ratio = cc4 / (cc4 - cc3);
    *point = Vec2d(_start.x() * ratio + _end.x() * (1.0 - ratio),
                   _start.y() * ratio + _end.y() * (1.0 - ratio));
    return true;
  }

  // 计算一个向量在线段上的内积
  double project_onto_unit(const Vec2d& point) const {
    return _unit_direction.inner_prod(point - _start);
  }

  // 计算一个向量在线段上的叉积
  double product_onto_unit(const Vec2d& point) const {
    return _unit_direction.cross_prod(point - _start);
  }

  // 计算直线上从线段扩展的二维点的垂直脚。
  double get_perpendicular_foot(const Vec2d& point, Vec2d* const foot_point) const {
    assert(foot_point);
    if (_length <= math_epsilon) {
      *foot_point = _start;
      return point.distance_to(_start);
    }
    const double x0 = point.x() - _start.x();
    const double y0 = point.y() - _start.y();
    const double proj = x0 * _unit_direction.x() + y0 * _unit_direction.y();
    *foot_point = _start + _unit_direction * proj;
    return std::abs(x0 * _unit_direction.y() - y0 * _unit_direction.x());
  }

#ifdef MYMATH_DBG
  std::string str() const {
    std::ostringstream ss;
    ss << "segment2d ( start = " << _start.str() << "  end = " << _end.str() << " )";
    return ss.str();
  }
#endif

  //
  // Static Function
  //
  static bool is_with_in(double val, double bound1, double bound2) {
    if (bound1 > bound2) {
      std::swap(bound1, bound2);
    }
    return val >= bound1 - math_epsilon && val <= bound2 + math_epsilon;
  }

private:
  Vec2d _start;
  Vec2d _end;
  Vec2d _unit_direction;
  double _heading = 0.0;
  double _length = 0.0;
};

}}

#endif
