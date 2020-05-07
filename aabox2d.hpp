#ifndef MYMATH_AABOX2D_HPP
#define MYMATH_AABOX2D_HPP

#include "mymath_config.h"

#include <cmath>
#include <vector>
#include <cassert>
#include <algorithm>

#ifdef MYMATH_DBG
#include <string>
#include <sstream>
#endif

#include "vec2d.hpp"
#include "math_utils.hpp"

namespace mypilot {
namespace mymath {

// 在二维中实现（无向）轴对齐的边界框。
class AABox2d {
public:
  // 创建一个轴对齐的盒子，其长度和宽度为零在原点。
  AABox2d() = default;

  /*
   * 创建一个轴对齐的盒子，其长度和宽度为零在原点。
   * 
   * center : 中间点坐标
   * length,width : 长(x轴),宽(y轴)
   */
  AABox2d(const Vec2d& center, const double length, const double width) : 
    _center(center),
    _length(length),
    _width(width),
    _half_length(length / 2.0),
    _half_width(width / 2.0) {
    assert(_length > -math_epsilon);
    assert(_width > -math_epsilon);
  }

  /*
   * 创建盒子使用角坐标
   * 
   * one_corner : 一个角坐标
   * opposite_corner : 对面的角坐标
   */
  AABox2d(const Vec2d& one_corner, const Vec2d& opposite_corner) : 
    AABox2d((one_corner + opposite_corner) / 2.0,
    std::abs(one_corner.x() - opposite_corner.x()),
    std::abs(one_corner.y() - opposite_corner.y())) {}

  // 创建盒子，使用所有的点坐标,这里的所有点都在这个盒子里。
  explicit AABox2d(const std::vector<Vec2d>& points)  {
    assert(!points.empty());
    double min_x = points[0].x();
    double max_x = points[0].x();
    double min_y = points[0].y();
    double max_y = points[0].y();
    // 遍历所有点坐标，找出最大最小值
    for (const auto& point : points) {
      min_x = std::min(min_x, point.x());
      max_x = std::max(max_x, point.x());
      min_y = std::min(min_y, point.y());
      max_y = std::max(max_y, point.y());
    }

    _center = {(min_x + max_x) / 2.0, (min_y + max_y) / 2.0};
    _length = max_x - min_x;
    _width = max_y - min_y;
    _half_length = _length / 2.0;
    _half_width = _width / 2.0;
  }

  const Vec2d& center() const { return _center; }
  double center_x() const { return _center.x(); }
  double center_y() const { return _center.y(); }
  double length() const { return _length; }
  double width() const { return _width; }
  double half_length() const { return _half_length; }
  double half_width() const { return _half_width; }
  double area() const { return _length * _width; }
  double min_x() const { return _center.x() - _half_length; }
  double max_x() const { return _center.x() + _half_length; }
  double min_y() const { return _center.y() - _half_width; }
  double max_y() const { return _center.y() + _half_width; }

  // 获取所有角坐标
  void get_all_corners(std::vector<Vec2d>* const corners) const {
    assert(corners);
    corners->clear();
    corners->reserve(4);
    corners->emplace_back(_center.x() + _half_length, _center.y() - _half_width);
    corners->emplace_back(_center.x() + _half_length, _center.y() + _half_width);
    corners->emplace_back(_center.x() - _half_length, _center.y() + _half_width);
    corners->emplace_back(_center.x() - _half_length, _center.y() - _half_width);
  }

  // 判断给定点坐标是否在盒子里
  bool is_point_in(const Vec2d& point) const {
    return std::abs(point.x() - _center.x()) <= _half_length + math_epsilon &&
      std::abs(point.y() - _center.y()) <= _half_width + math_epsilon;
  }

  // 给定点是否在盒子的边框上
  bool is_point_on_boundary(const Vec2d& point) const {
    const double dx = std::abs(point.x() - _center.x());
    const double dy = std::abs(point.y() - _center.y());
    return (std::abs(dx - _half_length) <= math_epsilon &&
        dy <= _half_width + math_epsilon) ||
      (std::abs(dy - _half_width) <= math_epsilon &&
        dx <= _half_length + math_epsilon);
  }

  // 确定点到盒子的距离
  double distance_to(const Vec2d& point) const {
    const double dx = std::abs(point.x() - _center.x()) - _half_length;
    const double dy = std::abs(point.y() - _center.y()) - _half_width;
    if (dx <= 0.0) {
      return std::max(0.0, dy);
    }
    if (dy <= 0.0) {
      return dx;
    }
    return hypot(dx, dy);
  }

  // 两个盒子的距离
  double distance_to(const AABox2d& box) const {
    const double dx =
      std::abs(box.center_x() - _center.x()) - box.half_length() - _half_length;
    const double dy =
      std::abs(box.center_y() - _center.y()) - box.half_width() - _half_width;
    if (dx <= 0.0) {
      return std::max(0.0, dy);
    }
    if (dy <= 0.0) {
      return dx;
    }
    return hypot(dx, dy);
  }

  // 检查两个盒子是否重叠
  bool has_overlap(const AABox2d& box) const {
    return std::abs(box.center_x() - _center.x()) <=
        box.half_length() + _half_length &&
      std::abs(box.center_y() - _center.y()) <=
        box.half_width() + _half_width;
  }

  // 通过给定偏移'shift_vec',移动盒子的中心。
  void shift(const Vec2d& shift_vec) { _center += shift_vec; }

  // 合并两个盒子
  void merge_from(const AABox2d& other_box) {
    const double x1 = std::min(min_x(), other_box.min_x());
    const double x2 = std::max(max_x(), other_box.max_x());
    const double y1 = std::min(min_y(), other_box.min_y());
    const double y2 = std::max(max_y(), other_box.max_y());
    _center = Vec2d((x1 + x2) / 2.0, (y1 + y2) / 2.0);
    _length = x2 - x1;
    _width = y2 - y1;
    _half_length = _length / 2.0;
    _half_width = _width / 2.0;
  }

  // 强制一个盒子包含一个点,扩大盒子。
  void merge_from(const Vec2d& other_point) {
    const double x1 = std::min(min_x(), other_point.x());
    const double x2 = std::max(max_x(), other_point.x());
    const double y1 = std::min(min_y(), other_point.y());
    const double y2 = std::max(max_y(), other_point.y());
    _center = Vec2d((x1 + x2) / 2.0, (y1 + y2) / 2.0);
    _length = x2 - x1;
    _width = y2 - y1;
    _half_length = _length / 2.0;
    _half_width = _width / 2.0;
  }

#if MYMATH_DBG
  std::string str() const {
    std::ostringstream ss;
    ss << "aabox2d ( center = " << _center.str() 
       << "  length = " << _length << "  width = " << _width
       << " )";
    return ss.str();
  }
#endif

private:
  Vec2d _center;
  double _length = 0.0;
  double _width = 0.0;
  double _half_length = 0.0;
  double _half_width = 0.0;
};

}}

#endif
