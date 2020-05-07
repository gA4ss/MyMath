#ifndef MYMATH_EULER_ANGLES_ZXY_HPP
#define MYMATH_EULER_ANGLES_ZXY_HPP

#include <cmath>
#include <Eigen/Geometry>

#include "math_utils.hpp"

namespace mypilot {
namespace mymath {

/*
 * 通过围绕正交坐标系的轴组成三个旋转，可以实现物体在3-D空间上的任何定向。 
 * 如果假定轴是不动的，则这些旋转被称为是外在的，否则为固有的。在这里，我们使用固有参照
 * 它相对于汽车的方向。 
 * 我们的车辆参考系遵循NovAtel的约定：
 * 轴x / y / z分别为右/前/上（RFU）。 
 * 特别是，我们通过三个角度来描述汽车的方向：
 * 1）（-pi / 2，pi / 2）中的螺旋对应于围绕x轴的旋转； 
 * 2）[-pi，pi）中的滚动对应于绕y轴的旋转； 
 * 3）[-pi，pi）中的偏航对应于围绕z轴的旋转。 
 * 当汽车水平时，俯仰为零；当机头朝上时，俯仰为正。 
 * 当汽车水平时，侧倾为零；当左侧部分向上时，侧倾为正。 
 * 当汽车朝北时，偏航为零，而朝西时，偏航为正。 
 * 反过来，在世界框架中，x / y / z轴指向东/北/上（ENU）。 
 * 这些角度表示从世界到车架的旋转。
 * 
 * T模板接受 double 与 float类型
 */
template <typename T>
class EulerAnglesZXY {
public:
  EulerAnglesZXY() : _roll(0), _pitch(0), _yaw(0) {}
  explicit EulerAnglesZXY(T yaw) : _roll(0), _pitch(0), _yaw(yaw) {}
  EulerAnglesZXY(T roll, T pitch, T yaw) : _roll(roll), _pitch(pitch), _yaw(yaw) {}
  // 使用四元组构造一个旋转
  EulerAnglesZXY(T qw, T qx, T qy, T qz) : 
    _roll(std::atan2(2.0 * (qw * qy - qx * qz), 
          2.0 * (square<T>(qw) + square<T>(qz)) - 1.0)),
          _pitch(std::asin(2.0 * (qw * qx + qy * qz))),
          _yaw(std::atan2(2.0 * (qw * qz - qx * qy), 
          2.0 * (square<T>(qw) + square<T>(qy)) - 1.0)) {}
  explicit EulerAnglesZXY(const Eigen::Quaternion<T>& q) : 
    EulerAnglesZXY(q.w(), q.x(), q.y(), q.z()) {}

  T roll() const { return _roll; }
  T pitch() const { return _pitch; }
  T yaw() const { return _yaw; }

  // 规范化 _roll, _pitch, _yaw 到 [-PI, PI)。
  void normalize() {
    _roll = normalize_angle(_roll);
    _pitch = normalize_angle(_pitch);
    _yaw = normalize_angle(_yaw);
  }

  // 验证一个指定的旋转角度是否有效，检查'pitch'角度是否在[-2PI, 2PI]之间。
  bool is_valid() {
    normalize();
    return _pitch < M_PI_2 && _pitch > -M_PI_2;
  }

  // 转换到一个Eigen的四元组
  Eigen::Quaternion<T> to_quaternion() const {
    T r = _roll * 0.5;
    T p = _pitch * 0.5;
    T y = _yaw * 0.5;

    T sr = std::sin(r);
    T sp = std::sin(p);
    T sy = std::sin(y);

    T cr = std::cos(r);
    T cp = std::cos(p);
    T cy = std::cos(y);

    T qw = cr * cp * cy - sr * sp * sy;
    T qx = cr * sp * cy - sr * cp * sy;
    T qy = cr * sp * sy + sr * cp * cy;
    T qz = cr * cp * sy + sr * sp * cy;

    // 方向相反
    if (qw < 0.0) return {-qw, -qx, -qy, -qz};
    return {qw, qx, qy, qz};
  }

private:
  T _roll;
  T _pitch;
  T _yaw;
};

using EulerAnglesZXYf = EulerAnglesZXY<float>;
using EulerAnglesZXYd = EulerAnglesZXY<double>;

}}

#endif
