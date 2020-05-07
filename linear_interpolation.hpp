#ifndef MYMATH_LINEAR_INTERPOLATION_HPP
#define MYMATH_LINEAR_INTERPOLATION_HPP

#include <cmath>

#include "math_utils.hpp"

#ifdef USE_PROTOC
#include <pnc_point.pb.h>
#else
#include "my_sl_point.hpp"
#include "my_trajectory_point.hpp"
#endif

namespace mypilot {
namespace mymath {

/*
 * 在两点之间进行线性插值。
 *
 * x0 : 第一点的坐标。
 * t0 : 第一点的插值参数。 
 * x1 : 第二点的坐标。
 * t1 : 第二点的插值参数。
 * t : 插值的插值参数。 
 * 返回插值点的坐标。
 */
template <typename T>
T lerp(const T& x0, const double t0, const T& x1, const double t1,
       const double t) {
  if (std::abs(t1 - t0) <= 1.0e-6) {
    // input time difference is too small
    return x0;
  }
  const double r = (t - t0) / (t1 - t0);
  const T x = x0 + r * (x1 - x0);
  return x;
}

/*
 * 两个角度之间的球面线性插值。这两个角度在[-M_PI，M_PI）范围内。
 * a0 : 第一个角度的值。
 * t0 : 第一个角度的插补参数。
 * a1 : 第二角度的值。
 * t1 : 第二角度的插补参数。
 * t : 用于插补的插补参数。
 * 返回球面插补角的值。
 */
double slerp(const double a0, const double t0, const double a1, const double t1,
             const double t) {
  if (std::abs(t1 - t0) <= math_epsilon) {
    // input time difference is too small
    return normalize_angle(a0);
  }
  const double a0_n = normalize_angle(a0);
  const double a1_n = normalize_angle(a1);
  double d = a1_n - a0_n;
  if (d > M_PI) {
    d = d - 2 * M_PI;
  } else if (d < -M_PI) {
    d = d + 2 * M_PI;
  }

  const double r = (t - t0) / (t1 - t0);
  const double a = a0_n + d * r;
  return normalize_angle(a);
}

SLPoint interpolate_using_linear_approximation(const SLPoint& p0,
                                               const SLPoint& p1, const double w) {
  assert(w >= 0.0);

  SLPoint p;
  p.set_s((1 - w) * p0.s() + w * p1.s());
  p.set_l((1 - w) * p0.l() + w * p1.l());
  return p;
}

PathPoint interpolate_using_linear_approximation(const PathPoint& p0,
                                                 const PathPoint& p1,
                                                 const double s) {
  double s0 = p0.s();
  double s1 = p1.s();
  assert(s0 <= s1);

  PathPoint path_point;
  double weight = (s - s0) / (s1 - s0);
  double x = (1 - weight) * p0.x() + weight * p1.x();
  double y = (1 - weight) * p0.y() + weight * p1.y();
  double theta = slerp(p0.theta(), p0.s(), p1.theta(), p1.s(), s);
  double kappa = (1 - weight) * p0.kappa() + weight * p1.kappa();
  double dkappa = (1 - weight) * p0.dkappa() + weight * p1.dkappa();
  double ddkappa = (1 - weight) * p0.ddkappa() + weight * p1.ddkappa();
  path_point.set_x(x);
  path_point.set_y(y);
  path_point.set_theta(theta);
  path_point.set_kappa(kappa);
  path_point.set_dkappa(dkappa);
  path_point.set_ddkappa(ddkappa);
  path_point.set_s(s);
  return path_point;
}

TrajectoryPoint interpolate_using_linear_approximation(const TrajectoryPoint& tp0,
                                                       const TrajectoryPoint& tp1,
                                                       const double t) {
  if (!tp0.has_path_point() || !tp1.has_path_point()) {
    TrajectoryPoint p;
#ifdef USE_PROTOC
    p.mutable_path_point()->CopyFrom(PathPoint());
#else
    p.set_path_point(PathPoint());
#endif
    return p;
  }
  const PathPoint pp0 = tp0.path_point();
  const PathPoint pp1 = tp1.path_point();
  double t0 = tp0.relative_time();
  double t1 = tp1.relative_time();

  TrajectoryPoint tp;
  tp.set_v(lerp(tp0.v(), t0, tp1.v(), t1, t));
  tp.set_a(lerp(tp0.a(), t0, tp1.a(), t1, t));
  tp.set_relative_time(t);

  PathPoint* path_point = tp.mutable_path_point();
  path_point->set_x(lerp(pp0.x(), t0, pp1.x(), t1, t));
  path_point->set_y(lerp(pp0.y(), t0, pp1.y(), t1, t));
  path_point->set_theta(slerp(pp0.theta(), t0, pp1.theta(), t1, t));
  path_point->set_kappa(lerp(pp0.kappa(), t0, pp1.kappa(), t1, t));
  path_point->set_dkappa(lerp(pp0.dkappa(), t0, pp1.dkappa(), t1, t));
  path_point->set_ddkappa(lerp(pp0.ddkappa(), t0, pp1.ddkappa(), t1, t));
  path_point->set_s(lerp(pp0.s(), t0, pp1.s(), t1, t));

  return tp;
}

}}

#endif
