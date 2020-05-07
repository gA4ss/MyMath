#ifndef MYMATH_VEC2D_HPP
#define MYMATH_VEC2D_HPP

#include "mymath_config.h"

#include <cassert>
#include <cmath>

#ifdef MYMATH_DBG
#include <string>
#include <sstream>
#endif

namespace mypilot {
namespace mymath {

static constexpr double math_epsilon = 1e-10;

class Vec2d {
public:
  constexpr Vec2d(const double x, const double y) noexcept : _x(x), _y(y) {}
  constexpr Vec2d() noexcept : Vec2d(0, 0) {}

  // 依照给定角度创建单位点
  static Vec2d create_unit_vec2d(const double angle) {
    return Vec2d(cos(angle), sin(angle));
  }

  double x() const { return _x; }
  double y() const { return _y; }
  void set_x(const double x) { _x = x; }
  void set_y(const double y) { _y = y; }

  double length() const {
    // 计算\sqrt(x^2 + y^2)
    return std::hypot(_x, _y);
  }

  double length_square() const {
    return _x * _x + _y * _y;
  }

  // 获取向量与x正半轴之间的夹角(返回弧度)
  double angle() const {
    return std::atan2(_y, _x);
  }

  void normalize() {
    const double l = length();
    if (l > math_epsilon) {
      _x /= l;
      _y /= l;
    }
  }

  // 计算与点other的距离
  double distance_to(const Vec2d& other) const {
    return std::hypot(_x - other._x, _y - other._y);
  }

  // 计算与点other的距离的平方
  double distance_square_to(const Vec2d& other) const {
    const double dx = _x - other._x;
    const double dy = _y - other._y;
    return dx * dx + dy * dy;
  }

  // 计算两个向量的叉积
  double cross_prod(const Vec2d& other) const {
    return _x * other.y() - _y * other.x();
  }

  // 求两点的内积
  double inner_prod(const Vec2d& other) const {
    return _x * other.x() + _y * other.y();
  }

  // 返回旋转angle度后的点
  Vec2d rotate(const double angle) const {
    return Vec2d(_x * cos(angle) - _y * sin(angle),
                 _x * sin(angle) + _y * cos(angle));
  }

  Vec2d operator+(const Vec2d& other) const {
    return Vec2d(_x + other.x(), _y + other.y());
  }

  Vec2d operator-(const Vec2d& other) const {
    return Vec2d(_x - other.x(), _y - other.y());
  }

  Vec2d operator*(const double ratio) const {
    return Vec2d(_x * ratio, _y * ratio);
  }

  Vec2d operator/(const double ratio) const {
    assert(std::abs(ratio) > math_epsilon);
    return Vec2d(_x / ratio, _y / ratio);
  }

  Vec2d& operator+=(const Vec2d& other) {
    _x += other.x();
    _y += other.y();
    return *this;
  }

  Vec2d& operator-=(const Vec2d& other) {
    _x -= other.x();
    _y -= other.y();
    return *this;
  }

  Vec2d& operator*=(const double ratio) {
    _x *= ratio;
    _y *= ratio;
    return *this;
  }

  Vec2d& operator/=(const double ratio) {
    assert(std::abs(ratio) > math_epsilon);
    _x /= ratio;
    _y /= ratio;
    return *this;
  }

  bool operator==(const Vec2d& other) const {
    return (std::abs(_x - other.x()) < math_epsilon &&
            std::abs(_y - other.y()) < math_epsilon);
  }

#ifdef MYMATH_DBG
  std::string str() const {
    std::ostringstream ss;
    ss << "vec2d ( x = " << _x << "  y = " << _y << " )";
    return ss.str();
  }
#endif

 protected:
  double _x = 0.0;
  double _y = 0.0;
};

// 一个点乘以一个标量
Vec2d operator*(const double ratio, const Vec2d& vec) {
  return vec * ratio;
}

}}

#endif
