#ifndef MYMATH_ANGLE_HPP
#define MYMATH_ANGLE_HPP

#include "mymath_config.h"

#include <cmath>
#include <cstdint>
#include <limits>

#ifdef USE_SIN_TABLE
#include "sin_table.h"
#endif

namespace mypilot {
namespace mymath {

/*
 * Angle类使用整数表示角度，并支持常用的操作，例如加法和减法，以及三角函数的使用。 
 * 具有专门的Angle类可防止代码重复，即用于某些任务，
 * 例如:计算角度差或以某个指定间隔（通常为[-pi，pi]）查找等效角度。以整数表示角度具有以下优点：
 * 1）精度控制的级别更高（'<' 表示 "精度低于"）：
 * Angle8 < Angle16 < float < Angle32 < double < Angle64 
 * 2）Angle8和Angle16允许通过64-KiB查找表实现超快速三角函数。 
 * 3）在相同的显示尺寸下精度更高。 
 * 鼓励使用Angle类。特别是应将Angle32用于纬度/经度（误差小于1cm）。
 * Angle16足够精确，可以进行定位/物体检测。
 * 
 * T模板接受 signed integer
 */
template <typename T>
class Angle {
public:
  static_assert(std::numeric_limits<T>::is_integer &&
                std::numeric_limits<T>::is_signed,
                "T must be a signed integer type");

  // 构造一个Angle对象从一个角度值
  static Angle from_deg(double value) {
    Angle a(std::lround(value * DEG_TO_RAW));
    return a;
  }

  // 构造一个Angle对象从一个弧度值
  static Angle from_rad(double value) {
    Angle a(std::lround(value * RAD_TO_RAW));
    return a;
  }

  // 构造一个Angle对象从一个角度值
  static double raw_to_deg(double raw) {
    return raw * RAW_TO_DEG;
  }

  // 构造一个Angle对象从一个弧度值
  static double raw_to_rad(double raw) {
    return raw * RAW_TO_RAD;
  }

  explicit Angle(T value=0) : _value(value) {}

  static constexpr T RAW_PI = std::numeric_limits<T>::min();
  static constexpr T RAW_PI_2 = -(std::numeric_limits<T>::min() >> 1);
  static constexpr double DEG_TO_RAW = RAW_PI / -180.0;
  static constexpr double RAD_TO_RAW = RAW_PI * -M_1_PI;
  static constexpr double RAW_TO_DEG = -180.0 / RAW_PI;
  static constexpr double RAW_TO_RAD = -M_PI / RAW_PI;

  T raw() const { return _value; }
  double to_deg() const { return _value * RAW_TO_DEG; }
  double to_rad() const { return _value * RAW_TO_RAD; }

  Angle operator+=(Angle other) {
    _value += other._value;
    return *this;
  }

  Angle operator-=(Angle other) {
    _value -= other._value;
    return *this;
  }

  template <typename Scalar>
  Angle operator*=(Scalar s) {
    _value = std::lround(_value * s);
    return *this;
  }

  template <typename Scalar>
  Angle operator/=(Scalar s) {
    _value = std::lround(_value / s);
    return *this;
  }

private:
  T _value;
};

using Angle8 = Angle<std::int8_t>;
using Angle16 = Angle<std::int16_t>;
using Angle32 = Angle<std::int32_t>;
using Angle64 = Angle<std::int64_t>;

template <typename T>
Angle<T> operator+(Angle<T> lhs, Angle<T> rhs) {
  lhs += rhs;
  return lhs;
}

template <typename T>
Angle<T> operator-(Angle<T> lhs, Angle<T> rhs) {
  lhs -= rhs;
  return lhs;
}

template <typename T, typename Scalar>
Angle<T> operator*(Angle<T> lhs, Scalar rhs) {
  lhs *= rhs;
  return lhs;
}

template <typename T, typename Scalar>
Angle<T> operator*(Scalar lhs, Angle<T> rhs) {
  rhs *= lhs;
  return rhs;
}

template <typename T, typename Scalar>
Angle<T> operator/(Angle<T> lhs, Scalar rhs) {
  lhs /= rhs;
  return lhs;
}

template <typename T>
double operator/(Angle<T> lhs, Angle<T> rhs) {
  return static_cast<double>(lhs.raw()) / rhs.raw();
}

template <typename T>
bool operator==(Angle<T> lhs, Angle<T> rhs) {
  return lhs.raw() == rhs.raw();
}

template <typename T>
bool operator!=(Angle<T> lhs, Angle<T> rhs) {
  return !(lhs == rhs);
}

float sin(Angle16 a) {
  auto idx = a.raw();
  if (idx < -Angle16::RAW_PI_2) {
    idx += Angle16::RAW_PI;
#ifdef USE_SIN_TABLE
    return -SIN_TABLE[idx % SIN_TABLE_SIZE];
#else
    return -std::sin(Angle16::raw_to_rad(idx));
#endif
  }
  if (idx < 0) {
#ifdef USE_SIN_TABLE
    return -SIN_TABLE[(-idx) % SIN_TABLE_SIZE];
#else
    return -std::sin(Angle16::raw_to_rad(idx));
#endif
  }
  if (idx < Angle16::RAW_PI_2) {
#ifdef USE_SIN_TABLE
    return SIN_TABLE[idx % SIN_TABLE_SIZE];
#else
    return std::sin(Angle16::raw_to_rad(idx));
#endif
  }
  idx = Angle16::RAW_PI - idx;
#ifdef USE_SIN_TABLE
  return SIN_TABLE[idx % SIN_TABLE_SIZE];
#else
    return std::sin(Angle16::raw_to_rad(idx));
#endif
}

float cos(Angle16 a) {
  Angle16 b(Angle16::RAW_PI_2 - a.raw());
  return sin(b);
}

float tan(Angle16 a) { return sin(a) / cos(a); }

float sin(Angle8 a) {
  Angle16 b(a.raw() << 8);
  return sin(b);
}

float cos(Angle8 a) {
  Angle16 b(a.raw() << 8);
  return cos(b);
}

float tan(Angle8 a) {
  Angle16 b(a.raw() << 8);
  return tan(b);
}

}}

#endif
