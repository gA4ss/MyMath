#ifndef MYMATH_EXTENDED_KALMAN_FILTER_HPP
#define MYMATH_EXTENDED_KALMAN_FILTER_HPP

#include <functional>
#include <utility>

#include <Eigen/Dense>
#include "matrix_operations.hpp"

namespace mypilot {
namespace mymath {

template<typename T, unsigned int XN, unsigned int ZN, unsigned int UN>
class ExtendedKalmanFilter {
public:
  ExtendedKalmanFilter() {
    _F.setIdentity();
    _Q.setZero();
    _H.setIdentity();
    _R.setZero();

    _x.setZero();
    _P.setZero();
    _y.setZero();
    _S.setZero();
    _K.setZero();
  }

  ExtendedKalmanFilter(const Eigen::Matrix<T, XN, 1>& x,
    const Eigen::Matrix<T, XN, XN>& P) :
    ExtendedKalmanFilter() {
    set_state_estimate(x, P);
  }

  virtual ~ExtendedKalmanFilter() {}

  /* 
   * 设置初始状态信念分布。
   * 
   * x : 信念分布的期望。
   * P : 状态信念分布的协方差。
   */
  void set_state_estimate(const Eigen::Matrix<T, XN, 1>& x,
    const Eigen::Matrix<T, XN, XN>& P) {
    _x = x;
    _P = P;
  }

  // 在零控制下更改系统转换功能
  void set_transition_model(
    std::function<
      Eigen::Matrix<T, XN, 1>(const Eigen::Matrix<T, XN, 1>&,
      const Eigen::Matrix<T, UN, 1>&)> f,
      const Eigen::Matrix<T, XN, XN>& F) {
    _f = std::move(f);
    _F = F;
  }

  // 更改过渡噪声的协方差矩阵
  void set_transition_noise(const Eigen::Matrix<T, XN, XN>& Q) { _Q = Q; }

  // 更改观察矩阵，将状态映射到观察
  void set_observation_model(
    std::function<Eigen::Matrix<T, ZN, 1>(
      const Eigen::Matrix<T, XN, 1>&)> h,
      const Eigen::Matrix<T, ZN, XN> &H) {
    _h = std::move(h);
    _H = H;
  }

  // 更改观测噪声的协方差矩阵
  void set_observation_noise(const Eigen::Matrix<T, ZN, ZN>& R) { _R = R; }
  // 更改当前状态信念分布的协方差矩阵
  void set_state_covariance(const Eigen::Matrix<T, XN, XN> &P) { _P = P; }

  // 获取变换噪声的协方差矩阵
  const Eigen::Matrix<T, XN, XN>& get_transition_noise() const { return _Q; }
  // 获取当前信念分布的期望
  Eigen::Matrix<T, XN, 1> get_state_mean() const { return _x; }
  // 获取信念分布的协方差矩阵。
  Eigen::Matrix<T, XN, XN> GetStateCovariance() const { return _P; }

  // 给定控制输入'u'(默认为0)更新状态置信度分布。
  void Predict(const Eigen::Matrix<T, UN, 1>& u=Eigen::Matrix<T, UN, 1>::Zero()) {
    _x = _f(_x, u);
    _P = _F * _P * _F.transpose() + _Q;
  }

  // 给定观察值'z'更新状态置信分布。
  void correct(const Eigen::Matrix<T, ZN, 1>& z) {
    _y = z - _h(_x);
    _S = _H * _P * _H.transpose() + _R;
    _K = _P * _H.transpose() * pseudo_inverse<T, ZN>(_S);
    _x = _x + _K * _y;
    _P = (Eigen::Matrix<T, XN, XN>::Identity() - _K * _H) * _P;
  }

private:
  // 当前状态信念分布的期望
  Eigen::Matrix<T, XN, 1> _x;
  // 当前状态信念分布的协方差
  Eigen::Matrix<T, XN, XN> _P;

  // 变换函数
  std::function<Eigen::Matrix<T, XN, 1>(
    const Eigen::Matrix<T, XN, 1>&,
    const Eigen::Matrix<T, UN, 1>&)> _f;

  // 零控制下的状态变换矩阵
  Eigen::Matrix<T, XN, XN> _F;
  // 状态变换噪声矩阵的协方差
  Eigen::Matrix<T, XN, XN> _Q;
  // 观察者函数
  std::function<Eigen::Matrix<T, ZN, 1>(const Eigen::Matrix<T, XN, 1>&)> _h;
  // 观察值矩阵
  Eigen::Matrix<T, ZN, XN> _H;
  // 观察值矩阵的协方差
  Eigen::Matrix<T, ZN, ZN> _R;
  // 更新;标记为成员，以防止重新分配内存。
  Eigen::Matrix<T, ZN, 1> _y;
  // 更新协方差；标记为成员，以防止重新分配内存。
  Eigen::Matrix<T, ZN, ZN> _S;
  // 卡尔曼增益；标记为成员，以防止重新分配内存。
  Eigen::Matrix<T, XN, ZN> _K;
};

}}

#endif
