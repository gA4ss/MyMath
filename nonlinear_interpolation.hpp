#ifndef COMMON_MYMATH_NONLINEAR_INTERPOLATION_HPP
#define COMMON_MYMATH_NONLINEAR_INTERPOLATION_HPP

//#include "modules/common/proto/pnc_point.pb.h"

#include <cmath>
#include <cassert>

#include "hermite_spline.h"
#include "integral.hpp"
#include "math_utils.hpp"

namespace mypilot {
namespace mymath {

PathPoint spline_interpolate(const PathPoint& p0, 
                             const PathPoint& p1,
                             const double s) {
  double s0 = p0.s();
  double s1 = p1.s();
  assert(s0 <= s && s <= s1);

  double theta_diff = normalize_angle(p1.theta() - p0.theta());

  std::array<double, 3> gx0{{0.0, p0.kappa(), p0.dkappa()}};
  std::array<double, 3> gx1{{theta_diff, p1.kappa(), p1.dkappa()}};

  HermiteSpline<double, 5> geometry_spline(gx0, gx1, s0, s1);
  auto func_cos_theta = [&geometry_spline, &p0](const double s) {
    auto theta = geometry_spline.evaluate(0, s) + p0.theta();
    return std::cos(theta);
  };
  auto func_sin_theta = [&geometry_spline, &p0](const double s) {
    auto theta = geometry_spline.evaluate(0, s) + p0.theta();
    return std::sin(theta);
  };

  double x = p0.x() + integrate_by_gauss_legendre<5>(func_cos_theta, s0, s);
  double y = p0.y() + integrate_by_gauss_legendre<5>(func_sin_theta, s0, s);
  double theta = normalize_angle(geometry_spline.evaluate(0, s) + p0.theta());
  double kappa = geometry_spline.evaluate(1, s);
  double dkappa = geometry_spline.evaluate(2, s);
  double d2kappa = geometry_spline.evaluate(3, s);

  PathPoint p;
  p.set_x(x);
  p.set_y(y);
  p.set_theta(theta);
  p.set_kappa(kappa);
  p.set_dkappa(dkappa);
  p.set_ddkappa(d2kappa);
  p.set_s(s);
  return p;
}

TrajectoryPoint spline_interpolate(const TrajectoryPoint& tp0,
                                   const TrajectoryPoint& tp1,
                                   const double t) {
  if (std::fabs(tp1.path_point().s() - tp0.path_point().s()) < 1.0e-4) {
    return tp1;
  }

  const PathPoint& pp0 = tp0.path_point();
  const PathPoint& pp1 = tp1.path_point();
  double t0 = tp0.relative_time();
  double t1 = tp1.relative_time();

  std::array<double, 2> dx0{{tp0.v(), tp0.a()}};
  std::array<double, 2> dx1{{tp1.v(), tp1.a()}};
  HermiteSpline<double, 3> dynamic_spline(dx0, dx1, t0, t1);

  double s0 = 0.0;
  auto func_v = [&dynamic_spline](const double t) {
    return dynamic_spline.evaluate(0, t);
  };
  double s1 = integrate_by_gauss_legendre<5>(func_v, t0, t1);
  double s = integrate_by_gauss_legendre<5>(func_v, t0, t);

  if (std::fabs(tp0.path_point().s() - s1) < 1.0e-4) {
    return tp1;
  }

  double v = dynamic_spline.evaluate(0, t);
  double a = dynamic_spline.evaluate(1, t);

  std::array<double, 2> gx0{{pp0.theta(), pp0.kappa()}};
  std::array<double, 2> gx1{{pp1.theta(), pp1.kappa()}};
  HermiteSpline<double, 3> geometry_spline(gx0, gx1, s0, s1);
  auto func_cos_theta = [&geometry_spline](const double s) {
    auto theta = geometry_spline.evaluate(0, s);
    return std::cos(theta);
  };
  auto func_sin_theta = [&geometry_spline](const double s) {
    auto theta = geometry_spline.evaluate(0, s);
    return std::sin(theta);
  };

  double x = pp0.x() + integrate_by_gauss_legendre<5>(func_cos_theta, s0, s);
  double y = pp0.y() + integrate_by_gauss_legendre<5>(func_sin_theta, s0, s);
  double theta = geometry_spline.evaluate(0, s);
  double kappa = geometry_spline.evaluate(1, s);
  double dkappa = geometry_spline.evaluate(2, s);
  double d2kappa = geometry_spline.evaluate(3, s);

  TrajectoryPoint tp;
  tp.set_v(v);
  tp.set_a(a);
  tp.set_relative_time(t);

  PathPoint* path_point = tp.mutable_path_point();
  path_point->set_x(x);
  path_point->set_y(y);
  path_point->set_theta(theta);
  path_point->set_kappa(kappa);
  path_point->set_dkappa(dkappa);
  path_point->set_ddkappa(d2kappa);
  path_point->set_s(s);
  return tp;
}

}}

#endif
