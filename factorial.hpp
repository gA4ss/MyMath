#ifndef MYMATH_FACTORIAL_H
#define MYMATH_FACTORIAL_H

namespace mypilot {
namespace mymath {

template <int N>
struct Factorial {
  enum { value = N * Factorial<N - 1>::value };
};

template <>
struct Factorial<0> {
  enum { value = 1 };
};

}}

#endif
