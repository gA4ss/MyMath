#ifndef MYMATH_MY_TRAJECTORYPOINT_HPP
#define MYMATH_MY_TRAJECTORYPOINT_HPP

#include "my_path_point.hpp"

namespace mypilot {
namespace mymath {

class TrajectoryPoint {
public:
  TrajectoryPoint() : _has_path_point(false) {}
  virtual ~TrajectoryPoint(){}

  bool has_path_point() const { return _has_path_point; }

  double v() const { return _v; }
  double a() const { return _a; }
  double relative_time() const { return _relative_time; }
  const PathPoint& path_point() const { return _path_point; }
  PathPoint* mutable_path_point() { return &_path_point; }

  void set_v(double v) { _v = v; }
  void set_a(double a) { _a = a; }
  void set_relative_time(double t) { _relative_time = t; }
  void set_path_point(PathPoint p) { _path_point = p; _has_path_point = true; }

private:
  double _v;
  double _a;
  double _relative_time;
  PathPoint _path_point;
  bool _has_path_point;
};

}}

#endif