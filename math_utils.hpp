#ifndef MYMATH_MATH_UTILS_HPP
#define MYMATH_MATH_UTILS_HPP

#include <cmath>
#include <limits>
#include <utility>

#include "vec2d.hpp"

namespace mypilot {
namespace mymath {

double sqr(const double x) { return x * x; }

/*
 * 计算两个2D向量(同一个起始点，不同的结束点)的叉积
 * start_point : 起始点
 * end_point_1 : 结束点1
 * end_point_2 : 结束点2
 */
double cross_prod(const Vec2d& start_point, 
                  const Vec2d& end_point_1,
                  const Vec2d& end_point_2) {
  return (end_point_1 - start_point).cross_prod(end_point_2 - start_point);
}

/*
 * 计算两个2D向量(同一个起始点，不同的结束点)的内积
 * start_point : 起始点
 * end_point_1 : 结束点1
 * end_point_2 : 结束点2
 */
double inner_prod(const Vec2d& start_point, 
                  const Vec2d& end_point_1,
                  const Vec2d& end_point_2) {
  return (end_point_1 - start_point).inner_prod(end_point_2 - start_point);
}

// 计算两条向量的叉积，以原点(0,0)为起点
double cross_prod(const double x0, const double y0, 
                  const double x1, const double y1) {
  return x0 * y1 - x1 * y0;
}

// 计算两条向量的内积，以原点(0,0)为起点
double inner_prod(const double x0, const double y0, 
                  const double x1, const double y1) {
  return x0 * x1 + y0 * y1;
}

// 将角度缩放到2pi内的同样的角度内,如果角度为负则求出它相同的正值。
double wrap_angle(const double angle) {
  const double new_angle = std::fmod(angle, M_PI * 2.0);
  return new_angle < 0 ? new_angle + M_PI * 2.0 : new_angle;
}

// 将角度标准化到[-PI, PI)内。
double normalize_angle(const double angle) {
  double a = std::fmod(angle + M_PI, 2.0 * M_PI);
  if (a < 0.0) {
    a += (2.0 * M_PI);
  }
  return a - M_PI;
}

// 计算从from到to角度之间的差距，并正则化到[0, PI)之间。
double angle_diff(const double from, const double to) {
  return normalize_angle(to - from);
}

// 随机选出一个在[s,t]之间的整数。
int random_int(const int s, const int t, unsigned int rand_seed=1) {
  if (s >= t) return s;
  return s + rand_r(&rand_seed) % (t - s + 1);
}

double random_double(const double s, const double t, unsigned int rand_seed=1) {
  return s + (t - s) / 16383.0 * (rand_r(&rand_seed) & 16383);
}

double gaussian(const double u, const double std, const double x) {
  return (1.0 / std::sqrt(2 * M_PI * std * std)) *
    std::exp(-(x - u) * (x - u) / (2 * std * std));
}

double sigmoid(const double x) { return 1.0 / (1.0 + std::exp(-x)); }

template <typename T>
inline T square(const T value) { return value * value; }

// 如果value小于bound1则返回bound1,如果大于bound2则返回bound2,否则返回value
template <typename T>
T clamp(const T value, T bound1, T bound2) {
  if (bound1 > bound2) {
    std::swap(bound1, bound2);
  }

  if (value < bound1) {
    return bound1;
  } else if (value > bound2) {
    return bound2;
  }
  return value;
}

// 将轴1中的点（x0，y0）转换为轴2中的点（x1，y1），其中从轴1到轴2的角度为theta（逆时针）
void rotate_axis(const double theta, const double x0, const double y0,
                 double *x1, double *y1) {
  assert(x1);
  assert(y1);

  const double cos_theta = std::cos(theta);
  const double sin_theta = std::sin(theta);
  *x1 = x0 * cos_theta + y0 * sin_theta;
  *y1 = -x0 * sin_theta + y0 * cos_theta;
}

inline std::pair<double, double> rfu_to_flu(const double x, const double y) {
  return std::make_pair(y, -x);
}

inline std::pair<double, double> flu_to_rfu(const double x, const double y) {
  return std::make_pair(-y, x);
}

// 求feat_dim维的向量feat_data进行l2范数
inline void l2_norm(int feat_dim, float *feat_data) {
  if (feat_dim == 0) {
    return;
  }
  // 正则化
  float l2norm = 0.0;
  for (int i = 0; i < feat_dim; ++i) {
    l2norm += feat_data[i] * feat_data[i];
  }
  if (l2norm == 0) {
    float val = 1.0 / std::sqrt(feat_dim);
    for (int i = 0; i < feat_dim; ++i) {
      feat_data[i] = val;
    }
  } else {
    l2norm = std::sqrt(l2norm);
    for (int i = 0; i < feat_dim; ++i) {
      feat_data[i] /= l2norm;
    }
  }
}

}}

#endif
