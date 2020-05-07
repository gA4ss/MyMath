#ifndef MYMATH_QUATERNION_HPP
#define MYMATH_QUATERNION_HPP

#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "euler_angles_zxy.hpp"
#include "math_utils.hpp"

#ifdef USE_PROTOC
#include <geometry.pb.h>
#else
#include "my_quaternion.hpp"
#endif

namespace mypilot {
namespace mymath {

/*
 * 返回[-PI, PI)区间的航向角度(弧度表示)
 * x/y/z 的方向分别是 东/北/上(左手坐标系)
 * 
 * qw : 四元组的'w'坐标
 * qx : 四元组的'x'坐标
 * qy : 四元组的'y'坐标
 * qz : 四元组的'z'坐标
 */
inline double quaternion_to_heading(const double qw, const double qx,
                                    const double qy, const double qz) {
  EulerAnglesZXYd euler_angles(qw, qx, qy, qz);
  // 当汽车指向北时,euler_angles.yaw()为0。但当汽车指向东时，航向为零
  return normalize_angle(euler_angles.yaw() + M_PI_2);
}

// 返回[-PI，PI）中的航向角（以弧度为单位），0为指向东方。请注意:"x/y/z为东/北/上"。
template <typename T>
inline T quaternion_to_heading(const Eigen::Quaternion<T> &q) {
  return quaternion_to_heading(q.w(), q.x(), q.y(), q.z());
}

/*
 * 返回一个四元数，该四元数以零滚动，零俯仰和指定的航向/偏航编码旋转。 
 * 请注意，航向在东方为零，偏航在北为零。 
 * 满足`quaternion_to_heading(heading_to_quaternion(h)) = h`
*/
template <typename T>
inline Eigen::Quaternion<T> heading_to_quaternion(T heading) {
  // 航向在东方为零，偏航在北为零
  EulerAnglesZXY<T> euler_angles(heading - M_PI_2);
  return euler_angles.to_quaternion();
}

/*
 * 将四元数定义的旋转应用于给定向量。请注意:"x/y/z"为东/北/上。
 * 
 * orientation : 定向四元组
 * original : 原始向量(左手坐标系)
 */
inline Eigen::Vector3d quaternion_rotate(const Quaternion& orientation,
                                         const Eigen::Vector3d& original) {
  Eigen::Quaternion<double> quaternion(orientation.qw(), orientation.qx(),
                                       orientation.qy(), orientation.qz());
  return quaternion.toRotationMatrix() * original;
}

// 按quaternion_rotate定义旋转角度进行反向操作
inline Eigen::Vector3d inverse_quaternion_rotate(const Quaternion& orientation,
                                                 const Eigen::Vector3d& rotated) {
  Eigen::Quaternion<double> quaternion(orientation.qw(), orientation.qx(),
                                       orientation.qy(), orientation.qz());
  return quaternion.toRotationMatrix().inverse() * rotated;
}

}}

#endif
