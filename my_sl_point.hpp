#ifndef MYMATH_MY_SLPOINT_HPP
#define MYMATH_MY_SLPOINT_HPP

namespace mypilot {
namespace mymath {

class SLPoint {
public:
  SLPoint(){}
  virtual ~SLPoint(){}

  double s() const { return _s; }
  double l() const { return _l; }

  void set_s(double s) { _s = s; }
  void set_l(double l) { _l = l; }

private:
  double _s;
  double _l;
};

}}

#endif