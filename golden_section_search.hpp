#ifndef MYMATH_GOLDEN_SECTION_SEARCH_HPP
#define MYMATH_GOLDEN_SECTION_SEARCH_HPP

#include <functional>

namespace mypilot {
namespace mymath {

/*
 * 给定一个在一定区间的一个单峰函数，寻找这个区间的最小值。
 * 参考 : https://en.wikipedia.org/wiki/Golden-section_search
 * func : 寻找的最小值的函数
 * lower_bound : 区间的下界
 * upper_bound : 区间的上界
 * tol : 误差的公差
 */
double golden_section_search(const std::function<double(double)> &func,
                             const double lower_bound, const double upper_bound,
                             const double tol = 1e-6) {
  constexpr double gr = 1.618033989;  // (sqrt(5) + 1) / 2

  double a = lower_bound;
  double b = upper_bound;

  double t = (b - a) / gr;      // 通过黄金分割比例点进行寻找
  double c = b - t;
  double d = a + t;

  // 循环遍历这个区间，直到对比出最小值
  while (std::abs(c - d) > tol) {
    if (func(c) < func(d)) {
      b = d;
    } else {
      a = c;
    }
    // 重新计算区域
    t = (b - a) / gr;
    c = b - t;
    d = a + t;
  }
  // 返回找到区域的中间点
  return (a + b) * 0.5;
}

}}

#endif
