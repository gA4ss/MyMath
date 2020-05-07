#ifndef MYMATH_BOX2D_HPP
#define MYMATH_BOX2D_HPP

#include <cmath>
#include <limits>
#include <vector>
#include <cassert>
#include <utility>
#include <algorithm>

#ifdef MYMATH_DBG
#include <string>
#include <sstream>
#endif

#include "math_utils.hpp"
#include "vec2d.hpp"
#include "line_segment2d.hpp"
#include "aabox2d.hpp"

namespace mypilot {
namespace mymath {

/* 
 * 二维矩形（无向）边界框
 * 尽管我们在该项目中使用偏航的约定（在East永久设置为0）
 * 迫使我们假定此处的'X/Y'，盒子为'东/北'（左手坐标系）。
 * 为了消除歧义，我们将与方向平行的矩形的轴称为'航向轴'。
 * 航向轴的大小称为'长度'，而垂直于该轴的轴的大小称为'宽度'。
 */
class Box2d {
public:
  Box2d() = default;

  // heading x轴与航向轴的夹角, length 航向轴长度, width 轴的垂足到航向轴的长度。
  Box2d(const Vec2d& center, const double heading, const double length,
        const double width) : 
    _center(center),
    _length(length),
    _width(width),
    _half_length(length / 2.0),
    _half_width(width / 2.0),
    _heading(heading),
    _cos_heading(cos(heading)),
    _sin_heading(sin(heading)) {
    assert(_length > -math_epsilon);
    assert(_width > -math_epsilon);
    init_corners();
  }

  // 盒子的宽度，垂直于航向轴方向。
  Box2d(const LineSegment2d& axis, const double width) : 
    _center(axis.center()),
    _length(axis.length()),
    _width(width),
    _half_length(axis.length() / 2.0),
    _half_width(width / 2.0),
    _heading(axis.heading()),
    _cos_heading(axis.cos_heading()),
    _sin_heading(axis.sin_heading()) {
    assert(_length > -math_epsilon);
    assert(_width > -math_epsilon);
    init_corners();
  }

  void init_corners() {
    const double dx1 = _cos_heading * _half_length;
    const double dy1 = _sin_heading * _half_length;
    const double dx2 = _sin_heading * _half_width;
    const double dy2 = -_cos_heading * _half_width;
    _corners.clear();
    _corners.emplace_back(_center.x() + dx1 + dx2, _center.y() + dy1 + dy2);
    _corners.emplace_back(_center.x() + dx1 - dx2, _center.y() + dy1 - dy2);
    _corners.emplace_back(_center.x() - dx1 - dx2, _center.y() - dy1 - dy2);
    _corners.emplace_back(_center.x() - dx1 + dx2, _center.y() - dy1 + dy2);

    for (auto& corner : _corners) {
      _max_x = std::fmax(corner.x(), _max_x);
      _min_x = std::fmin(corner.x(), _min_x);
      _max_y = std::fmax(corner.y(), _max_y);
      _min_y = std::fmin(corner.y(), _min_y);
    }
  }

  // 通过轴对齐盒子构造
  explicit Box2d(const AABox2d& aabox) : 
    _center(aabox.center()),
    _length(aabox.length()),
    _width(aabox.width()),
    _half_length(aabox.half_length()),
    _half_width(aabox.half_width()),
    _heading(0.0),
    _cos_heading(1.0),
    _sin_heading(0.0) {
    assert(_length > -math_epsilon);
    assert(_width > -math_epsilon);
  }

  // 使用两个对应的顶点坐标创建盒子
  static Box2d create_aabox(const Vec2d& one_corner,
                            const Vec2d& opposite_corner) {
    const double x1 = std::min(one_corner.x(), opposite_corner.x());
    const double x2 = std::max(one_corner.x(), opposite_corner.x());
    const double y1 = std::min(one_corner.y(), opposite_corner.y());
    const double y2 = std::max(one_corner.y(), opposite_corner.y());
    return Box2d({(x1 + x2) / 2.0, (y1 + y2) / 2.0}, 0.0, x2 - x1, y2 - y1);
  }

  const Vec2d& center() const { return _center; }
  double center_x() const { return _center.x(); }
  double center_y() const { return _center.y(); }
  double length() const { return _length; }
  double width() const { return _width; }
  double half_length() const { return _half_length; }
  double half_width() const { return _half_width; }
  double heading() const { return _heading; }
  double cos_heading() const { return _cos_heading; }
  double sin_heading() const { return _sin_heading; }
  double area() const { return _length * _width; }
  double diagonal() const { return std::hypot(_length, _width); }

  // 获取所有的顶点
  void get_all_corners(std::vector<Vec2d>* const corners) const {
    if (corners == nullptr) return;
    *corners = _corners;
  }
  std::vector<Vec2d> get_all_corners() const { return _corners; }

  // 按照'ccw'顺序对端点进行排序
  static void sort_corners_by_ccw(std::vector<Vec2d>& corners) {
    size_t num_points = corners.size();
    assert(num_points >= 3);

    // 确认点的顺序为'ccw'顺序
    double area = 0.0;
    for (int i = 1; i < num_points; ++i) {
      area += cross_prod(corners[0], corners[i - 1], corners[i]);
    }
    if (area < 0) {
      area = -area;
      std::reverse(corners.begin(), corners.end());
    }
    area /= 2.0;
    assert(area > math_epsilon);
    return;
  }

  std::vector<LineSegment2d> get_all_segments() const {
    std::vector<Vec2d> points = get_all_corners();
    return get_all_segments(points);
  }

  static std::vector<LineSegment2d> get_all_segments(std::vector<Vec2d>& points) {
    size_t num_points = points.size();
    assert(num_points >= 3);

    // 对端点进行排序
    sort_corners_by_ccw(points);

    // 构造线段
    std::vector<LineSegment2d> line_segments;
    line_segments.reserve(num_points);
    for (int i = 0; i < num_points; ++i) {
      int next_i = i >= num_points - 1 ? 0 : i + 1;
      line_segments.emplace_back(points[i], points[next_i]);
    }
    return line_segments;
  }

  // 判断端点是否在盒子内
  bool is_point_in(const Vec2d& point) const {
    const double x0 = point.x() - _center.x();
    const double y0 = point.y() - _center.y();
    const double dx = std::abs(x0 * _cos_heading + y0 * _sin_heading);
    const double dy = std::abs(-x0 * _sin_heading + y0 * _cos_heading);
    return dx <= _half_length + math_epsilon && dy <= _half_width + math_epsilon;
  }

  // 在点上在盒子的边界上
  bool is_point_on_boundary(const Vec2d& point) const {
    const double x0 = point.x() - _center.x();
    const double y0 = point.y() - _center.y();
    const double dx = std::abs(x0 * _cos_heading + y0 * _sin_heading);
    const double dy = std::abs(x0 * _sin_heading - y0 * _cos_heading);
    return (std::abs(dx - _half_length) <= math_epsilon &&
            dy <= _half_width + math_epsilon) ||
          (std::abs(dy - _half_width) <= math_epsilon &&
            dx <= _half_length + math_epsilon);
  }

  // 测试与点的距离
  double distance_to(const Vec2d& point) const {
    const double x0 = point.x() - _center.x();
    const double y0 = point.y() - _center.y();
    const double dx =
      std::abs(x0 * _cos_heading + y0 * _sin_heading) - _half_length;
    const double dy =
      std::abs(x0 * _sin_heading - y0 * _cos_heading) - _half_width;
    if (dx <= 0.0)
      return std::max(0.0, dy);
    if (dy <= 0.0)
      return dx;
    return hypot(dx, dy);
  }

  // 与线段的距离
  double distance_to(const LineSegment2d& line_segment) const {
    if (line_segment.length() <= math_epsilon)
      return distance_to(line_segment.start());

    const double ref_x1 = line_segment.start().x() - _center.x();
    const double ref_y1 = line_segment.start().y() - _center.y();
    double x1 = ref_x1 * _cos_heading + ref_y1 * _sin_heading;
    double y1 = ref_x1 * _sin_heading - ref_y1 * _cos_heading;
    double box_x = _half_length;
    double box_y = _half_width;
    int gx1 = (x1 >= box_x ? 1 : (x1 <= -box_x ? -1 : 0));
    int gy1 = (y1 >= box_y ? 1 : (y1 <= -box_y ? -1 : 0));
    if (gx1 == 0 && gy1 == 0)
      return 0.0;
    
    const double ref_x2 = line_segment.end().x() - _center.x();
    const double ref_y2 = line_segment.end().y() - _center.y();
    double x2 = ref_x2 * _cos_heading + ref_y2 * _sin_heading;
    double y2 = ref_x2 * _sin_heading - ref_y2 * _cos_heading;
    int gx2 = (x2 >= box_x ? 1 : (x2 <= -box_x ? -1 : 0));
    int gy2 = (y2 >= box_y ? 1 : (y2 <= -box_y ? -1 : 0));
    if (gx2 == 0 && gy2 == 0)
      return 0.0;

    if (gx1 < 0 || (gx1 == 0 && gx2 < 0)) {
      x1 = -x1;
      gx1 = -gx1;
      x2 = -x2;
      gx2 = -gx2;
    }
    if (gy1 < 0 || (gy1 == 0 && gy2 < 0)) {
      y1 = -y1;
      gy1 = -gy1;
      y2 = -y2;
      gy2 = -gy2;
    }
    if (gx1 < gy1 || (gx1 == gy1 && gx2 < gy2)) {
      std::swap(x1, y1);
      std::swap(gx1, gy1);
      std::swap(x2, y2);
      std::swap(gx2, gy2);
      std::swap(box_x, box_y);
    }
    if (gx1 == 1 && gy1 == 1) {
      switch (gx2 * 3 + gy2) {
        case 4:
          return ptseg_distance(box_x, box_y, x1, y1, x2, y2,
                                line_segment.length());
        case 3:
          return (x1 > x2) ? (x2 - box_x) : 
            ptseg_distance(box_x, box_y, x1, y1, x2, y2,
                           line_segment.length());
        case 2:
          return (x1 > x2) ? ptseg_distance(box_x, -box_y, x1, y1, x2, y2,
                                            line_segment.length())
                           : ptseg_distance(box_x, box_y, x1, y1, x2, y2,
                                            line_segment.length());
        case -1:
          return cross_prod({x1, y1}, {x2, y2}, {box_x, -box_y}) >= 0.0
                    ? 0.0
                    : ptseg_distance(box_x, -box_y, x1, y1, x2, y2,
                                     line_segment.length());
        case -4:
          return cross_prod({x1, y1}, {x2, y2}, {box_x, -box_y}) <= 0.0
                    ? ptseg_distance(box_x, -box_y, x1, y1, x2, y2,
                                     line_segment.length())
                    : (cross_prod({x1, y1}, {x2, y2}, {-box_x, box_y}) <= 0.0
                            ? 0.0
                            : ptseg_distance(-box_x, box_y, x1, y1, x2, y2,
                                             line_segment.length()));
      }/* end switch */
    } else {
      switch (gx2 * 3 + gy2) {
        case 4:
          return (x1 < x2) ? (x1 - box_x)
                           : ptseg_distance(box_x, box_y, x1, y1, x2, y2,
                                            line_segment.length());
        case 3:
          return std::min(x1, x2) - box_x;
        case 1:
        case -2:
          return cross_prod({x1, y1}, {x2, y2}, {box_x, box_y}) <= 0.0
                    ? 0.0
                    : ptseg_distance(box_x, box_y, x1, y1, x2, y2,
                                     line_segment.length());
        case -3:
          return 0.0;
      }
    }

    // 到达这里就出错
    assert(0);
    // unimplemented state: gx1, gy1, gx2 gy2
    return 0.0;
  }

  // 计算两个盒子的距离,重叠返回0。
  double distance_to(const Box2d& box) const {
    //return Polygon2d(box).distance_to(*this);
    std::vector<Vec2d> corners = box.get_all_corners();
    assert(corners.size() >= 3);
    sort_corners_by_ccw(corners);

    std::vector<LineSegment2d> line_segments = get_all_segments(corners);
    double distance = std::numeric_limits<double>::infinity();
    for (int i = 0; i < corners.size(); ++i) {
      // 检测box的顶点是否在当前盒子内，在则返回。
      if (is_point_in(corners[i])) return 0.0;
      distance = std::min(distance, box.distance_to(line_segments[i]));
    }
    return distance;
  }

  // 检查与线段是否重叠
  bool has_overlap(const LineSegment2d& line_segment) const {
    if (line_segment.length() <= math_epsilon) {
      return is_point_in(line_segment.start());
    }
    if (std::fmax(line_segment.start().x(), line_segment.end().x()) < min_x() ||
        std::fmin(line_segment.start().x(), line_segment.end().x()) > max_x() ||
        std::fmax(line_segment.start().y(), line_segment.end().y()) < min_x() ||
        std::fmin(line_segment.start().y(), line_segment.end().y()) > max_x()) {
      return false;
    }
    return distance_to(line_segment) <= math_epsilon;
  }

  // 检查两个盒子是否重合
  bool has_overlap(const Box2d& box) const {
    if (box.max_x() < min_x() || box.min_x() > max_x() || box.max_y() < min_y() ||
        box.min_y() > max_y()) {
      return false;
    }

    const double shift_x = box.center_x() - _center.x();
    const double shift_y = box.center_y() - _center.y();

    const double dx1 = _cos_heading * _half_length;
    const double dy1 = _sin_heading * _half_length;
    const double dx2 = _sin_heading * _half_width;
    const double dy2 = -_cos_heading * _half_width;
    const double dx3 = box.cos_heading() * box.half_length();
    const double dy3 = box.sin_heading() * box.half_length();
    const double dx4 = box.sin_heading() * box.half_width();
    const double dy4 = -box.cos_heading() * box.half_width();

    return std::abs(shift_x * _cos_heading + shift_y * _sin_heading) <=
           std::abs(dx3 * _cos_heading + dy3 * _sin_heading) +
           std::abs(dx4 * _cos_heading + dy4 * _sin_heading) +
              _half_length &&
           std::abs(shift_x * _sin_heading - shift_y * _cos_heading) <=
           std::abs(dx3 * _sin_heading - dy3 * _cos_heading) +
           std::abs(dx4 * _sin_heading - dy4 * _cos_heading) +
              _half_width &&
           std::abs(shift_x * box.cos_heading() + shift_y * box.sin_heading()) <=
           std::abs(dx1 * box.cos_heading() + dy1 * box.sin_heading()) +
           std::abs(dx2 * box.cos_heading() + dy2 * box.sin_heading()) +
              box.half_length() &&
           std::abs(shift_x * box.sin_heading() - shift_y * box.cos_heading()) <=
           std::abs(dx1 * box.sin_heading() - dy1 * box.cos_heading()) +
           std::abs(dx2 * box.sin_heading() - dy2 * box.cos_heading()) +
              box.half_width();
  }

  // 获取一个最小的轴对齐盒子
  AABox2d get_aabox() const {
    const double dx1 = std::abs(_cos_heading * _half_length);
    const double dy1 = std::abs(_sin_heading * _half_length);
    const double dx2 = std::abs(_sin_heading * _half_width);
    const double dy2 = std::abs(_cos_heading * _half_width);
    return AABox2d(_center, (dx1 + dx2) * 2.0, (dy1 + dy2) * 2.0);
  }

  // 以中央点进行'rotate_angle'角度旋转
  void rotate_from_center(const double rotate_angle) {
    _heading = normalize_angle(_heading + rotate_angle);
    _cos_heading = std::cos(_heading);
    _sin_heading = std::sin(_heading);
    init_corners();
  }

  // 通过给定的偏移'shift_vec'移动盒子
  void shift(const Vec2d& shift_vec) {
    _center += shift_vec;
    init_corners();
  }

  // 纵向地扩展盒子'extension_length'长度
  void longitudinal_extend(const double extension_length) {
    _length += extension_length;
    _half_length += extension_length / 2.0;
    init_corners();
  }

  void lateral_extend(const double extension_length) {
    _width += extension_length;
    _half_width += extension_length / 2.0;
    init_corners();
  }

#if MYMATH_DBG
  std::string str() const {
    std::ostringstream ss;
    ss << "box2d ( center = " << _center.str()
       << "  heading = " << _heading
       << "  length = " << _length
       << "  width = " << _width
       << " )";
  return ss.str();
  }
#endif

  double max_x() const { return _max_x; }
  double min_x() const { return _min_x; }
  double max_y() const { return _max_y; }
  double min_y() const { return _min_y; }

  static double ptseg_distance(double query_x, double query_y, double start_x,
                               double start_y, double end_x, double end_y,
                               double length) {
    const double x0 = query_x - start_x;
    const double y0 = query_y - start_y;
    const double dx = end_x - start_x;
    const double dy = end_y - start_y;
    const double proj = x0 * dx + y0 * dy;
    if (proj <= 0.0) {
      return hypot(x0, y0);
    }
    if (proj >= length * length) {
      return hypot(x0 - dx, y0 - dy);
    }
    return std::abs(x0 * dy - y0 * dx) / length;
  }

private:
  Vec2d _center;
  double _length = 0.0;
  double _width = 0.0;
  double _half_length = 0.0;
  double _half_width = 0.0;
  double _heading = 0.0;
  double _cos_heading = 1.0;
  double _sin_heading = 0.0;

  std::vector<Vec2d> _corners;

  double _max_x = std::numeric_limits<double>::min();
  double _min_x = std::numeric_limits<double>::max();
  double _max_y = std::numeric_limits<double>::min();
  double _min_y = std::numeric_limits<double>::max();
};

}}

#endif
