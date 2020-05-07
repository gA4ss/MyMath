#ifndef MYMATH_MY_PATHPOINT_HPP
#define MYMATH_MY_PATHPOINT_HPP

namespace mypilot {
namespace mymath {

class PathPoint {
public:
  PathPoint(){}
  virtual ~PathPoint(){}

  double x() const { return _x; }
  double y() const { return _y; }
  double theta() const { return _theta; }
  double kappa() const { return _kappa; }
  double dkappa() const { return _dkappa; }
  double ddkappa() const { return _ddkappa; }
  double s() const { return _s; }

  void set_x(double x) { _x = x; }
  void set_y(double y) { _y = y; }
  void set_theta(double theta) { _theta = theta; }
  void set_kappa(double kappa) { _kappa = kappa; }
  void set_dkappa(double dkappa) { _dkappa = dkappa; }
  void set_ddkappa(double ddkappa) { _ddkappa = ddkappa; }
  void set_s(double s) { _s = s; }

private:
  double _x;
  double _y;
  double _theta;
  double _kappa;
  double _dkappa;
  double _ddkappa;
  double _s;
};

}}

#endif