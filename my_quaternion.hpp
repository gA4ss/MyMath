#ifndef MYMATH_MY_QUATERNION_HPP
#define MYMATH_MY_QUATERNION_HPP

namespace mypilot {
namespace mymath {

class Quaternion {
public:
  Quaternion(){}
  virtual ~Quaternion(){}

  double qw() const { return _qw; }
  double qx() const { return _qx; }
  double qy() const { return _qy; }
  double qz() const { return _qz; }

  void set_qw(double w) { _qw = w; }
  void set_qx(double x) { _qx = x; }
  void set_qy(double y) { _qy = y; }
  void set_qz(double z) { _qz = z; }

private:
  double _qw;
  double _qx;
  double _qy;
  double _qz;
};

}}

#endif