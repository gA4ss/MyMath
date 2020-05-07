#ifndef MYMATH_INTEGRAL_HPP
#define MYMATH_INTEGRAL_HPP

#include <array>
#include <cassert>
#include <functional>
#include <utility>
#include <vector>

namespace mypilot {
namespace mymath {

/* 
 * simpson积分公式
 * 在平面直角坐标系中，由任意三点(x1, y1), (x2, y2), (x3, y3) (x1<x2<x3，x2 = (x1 + x3)/2)。
 * 确定的抛物线y=f(x)在[x1, x3]的定积分为"辛普森公式"。
 */
double integrate_by_simpson(const std::vector<double>& funv_vec, const double dx,
                            const std::size_t nsteps) {
  assert(1 == nsteps & 1);
  double sum1 = 0.0;
  double sum2 = 0.0;
  for (std::size_t i = 1; i + 1 < nsteps; ++i) {
    if ((i & 1) != 0) {
      sum1 += funv_vec[i];
    } else {
      sum2 += funv_vec[i];
    }
  }
  return dx / 3.0 * (4.0 * sum1 + 2.0 * sum2 + funv_vec[0] + funv_vec[nsteps - 1]);
}

// 梯形数值积分
double integrate_by_trapezoidal(const std::vector<double>& funv_vec,
                                const double dx, const std::size_t nsteps) {
  double sum = 0;
  for (std::size_t i = 1; i + 1 < nsteps; ++i) {
    sum += funv_vec[i];
  }
  return dx * sum + 0.5 * dx * (funv_vec[0] + funv_vec[nsteps - 1]);
}

// 获取不同有序Gauss-Legendre，积分的分数和权重。当前支持阶2-10。其他输入阶将触发编译错误。 
template <std::size_t N>
std::pair<std::array<double, N>, std::array<double, N>>
get_gauss_legendre_points();

template <>
inline std::pair<std::array<double, 2>, std::array<double, 2>>
get_gauss_legendre_points<2>() {
  std::array<double, 2> x;
  x[0] = -5.77350269189625764507e-01;
  x[1] = 5.77350269189625764507e-01;

  std::array<double, 2> w;
  w[0] = 1.0;
  w[1] = 1.0;

  return std::make_pair(x, w);
}

template <>
inline std::pair<std::array<double, 3>, std::array<double, 3>>
get_gauss_legendre_points<3>() {
  std::array<double, 3> x;
  x[0] = 0.00000000000000000000e+00;
  x[1] = 7.74596669241483377010e-01;
  x[2] = -7.74596669241483377010e-01;

  std::array<double, 3> w;
  w[0] = 8.88888888888888888877e-01;
  w[1] = 5.55555555555555555562e-01;
  w[2] = 5.55555555555555555562e-01;

  return std::make_pair(x, w);
}

template <>
inline std::pair<std::array<double, 4>, std::array<double, 4>>
get_gauss_legendre_points<4>() {
  std::array<double, 4> x;
  x[0] = 3.39981043584856264792e-01;
  x[1] = -3.39981043584856264792e-01;
  x[2] = 8.61136311594052575248e-01;
  x[3] = -8.61136311594052575248e-01;

  std::array<double, 4> w;
  w[0] = 6.52145154862546142644e-01;
  w[1] = 6.52145154862546142644e-01;
  w[2] = 3.47854845137453857383e-01;
  w[3] = 3.47854845137453857383e-01;

  return std::make_pair(x, w);
}

template <>
inline std::pair<std::array<double, 5>, std::array<double, 5>>
get_gauss_legendre_points<5>() {
  std::array<double, 5> x;
  x[0] = 0.00000000000000000000e+00;
  x[1] = 5.38469310105683091018e-01;
  x[2] = -5.38469310105683091018e-01;
  x[3] = 9.06179845938663992811e-01;
  x[4] = -9.06179845938663992811e-01;

  std::array<double, 5> w;
  w[0] = 5.68888888888888888883e-01;
  w[1] = 4.78628670499366468030e-01;
  w[2] = 4.78628670499366468030e-01;
  w[3] = 2.36926885056189087515e-01;
  w[4] = 2.36926885056189087515e-01;

  return std::make_pair(x, w);
}

template <>
inline std::pair<std::array<double, 6>, std::array<double, 6>>
get_gauss_legendre_points<6>() {
  std::array<double, 6> x;
  x[0] = 6.61209386466264513688e-01;
  x[1] = -6.61209386466264513688e-01;
  x[2] = 2.38619186083196908630e-01;
  x[3] = -2.38619186083196908630e-01;
  x[4] = 9.32469514203152027832e-01;
  x[5] = -9.32469514203152027832e-01;

  std::array<double, 6> w;
  w[0] = 3.60761573048138607569e-01;
  w[1] = 3.60761573048138607569e-01;
  w[2] = 4.67913934572691047389e-01;
  w[3] = 4.67913934572691047389e-01;
  w[4] = 1.71324492379170345043e-01;
  w[5] = 1.71324492379170345043e-01;

  return std::make_pair(x, w);
}

template <>
inline std::pair<std::array<double, 7>, std::array<double, 7>>
get_gauss_legendre_points<7>() {
  std::array<double, 7> x;
  x[0] = 0.00000000000000000000e+00;
  x[1] = 4.05845151377397166917e-01;
  x[2] = -4.05845151377397166917e-01;
  x[3] = 7.41531185599394439864e-01;
  x[4] = -7.41531185599394439864e-01;
  x[5] = 9.49107912342758524541e-01;
  x[6] = -9.49107912342758524541e-01;

  std::array<double, 7> w;
  w[0] = 4.17959183673469387749e-01;
  w[1] = 3.81830050505118944961e-01;
  w[2] = 3.81830050505118944961e-01;
  w[3] = 2.79705391489276667890e-01;
  w[4] = 2.79705391489276667890e-01;
  w[5] = 1.29484966168869693274e-01;
  w[6] = 1.29484966168869693274e-01;

  return std::make_pair(x, w);
}

template <>
inline std::pair<std::array<double, 8>, std::array<double, 8>>
get_gauss_legendre_points<8>() {
  std::array<double, 8> x;
  x[0] = 1.83434642495649804936e-01;
  x[1] = -1.83434642495649804936e-01;
  x[2] = 5.25532409916328985830e-01;
  x[3] = -5.25532409916328985830e-01;
  x[4] = 7.96666477413626739567e-01;
  x[5] = -7.96666477413626739567e-01;
  x[6] = 9.60289856497536231661e-01;
  x[7] = -9.60289856497536231661e-01;

  std::array<double, 8> w;
  w[0] = 3.62683783378361982976e-01;
  w[1] = 3.62683783378361982976e-01;
  w[2] = 3.13706645877887287338e-01;
  w[3] = 3.13706645877887287338e-01;
  w[4] = 2.22381034453374470546e-01;
  w[5] = 2.22381034453374470546e-01;
  w[6] = 1.01228536290376259154e-01;
  w[7] = 1.01228536290376259154e-01;

  return std::make_pair(x, w);
}

template <>
inline std::pair<std::array<double, 9>, std::array<double, 9>>
get_gauss_legendre_points<9>() {
  std::array<double, 9> x;
  x[0] = 0.00000000000000000000e+00;
  x[1] = 8.36031107326635794313e-01;
  x[2] = -8.36031107326635794313e-01;
  x[3] = 9.68160239507626089810e-01;
  x[4] = -9.68160239507626089810e-01;
  x[5] = 3.24253423403808929042e-01;
  x[6] = -3.24253423403808929042e-01;
  x[7] = 6.13371432700590397285e-01;
  x[8] = -6.13371432700590397285e-01;

  std::array<double, 9> w;
  w[0] = 3.30239355001259763154e-01;
  w[1] = 1.80648160694857404059e-01;
  w[2] = 1.80648160694857404059e-01;
  w[3] = 8.12743883615744119737e-02;
  w[4] = 8.12743883615744119737e-02;
  w[5] = 3.12347077040002840057e-01;
  w[6] = 3.12347077040002840057e-01;
  w[7] = 2.60610696402935462313e-01;
  w[8] = 2.60610696402935462313e-01;

  return std::make_pair(x, w);
}

template <>
inline std::pair<std::array<double, 10>, std::array<double, 10>>
get_gauss_legendre_points<10>() {
  std::array<double, 10> x;
  x[0] = 1.48874338981631210881e-01;
  x[1] = -1.48874338981631210881e-01;
  x[2] = 4.33395394129247190794e-01;
  x[3] = -4.33395394129247190794e-01;
  x[4] = 6.79409568299024406207e-01;
  x[5] = -6.79409568299024406207e-01;
  x[6] = 8.65063366688984510759e-01;
  x[7] = -8.65063366688984510759e-01;
  x[8] = 9.73906528517171720066e-01;
  x[9] = -9.73906528517171720066e-01;

  std::array<double, 10> w;
  w[0] = 2.95524224714752870187e-01;
  w[1] = 2.95524224714752870187e-01;
  w[2] = 2.69266719309996355105e-01;
  w[3] = 2.69266719309996355105e-01;
  w[4] = 2.19086362515982044000e-01;
  w[5] = 2.19086362515982044000e-01;
  w[6] = 1.49451349150580593150e-01;
  w[7] = 1.49451349150580593150e-01;
  w[8] = 6.66713443086881375920e-02;
  w[9] = 6.66713443086881375920e-02;

  return std::make_pair(x, w);
}

/*
 * 通过5阶Gauss-Legendre方法计算目标单变量函数的积分。
 * 从下限到上限，求给定目标函数以及上下限的积分，
 * 使用5阶Gauss-Legendre计算逼近积分，目标函数必须是平滑函数。
 * 
 * 例子:
 * 目标函数 : auto func = [](const double x) {return x * x;};
 *           double integral = gauss_legendre(func, -2, 3);
 * 这将返回在定义域[-2,3]的x^2函数的近似积分。
 */
template <std::size_t N>
double integrate_by_gauss_legendre(const std::function<double(double)>& func,
                                   const double lower_bound,
                                   const double upper_bound) {
  auto p = get_gauss_legendre_points<N>();

  std::array<double, N> x = p.first;
  std::array<double, N> w = p.second;

  const double t = (upper_bound - lower_bound) * 0.5;
  const double m = (upper_bound + lower_bound) * 0.5;

  double integral = 0.0;
  for (size_t i = 0; i < N; ++i) {
    integral += w[i] * func(t * x[i] + m);
  }

  return integral * t;
}

}}

#endif
