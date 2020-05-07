#ifndef MYMATH_KALMAN_FILTER_HPP
#define MYMATH_KALMAN_FILTER_HPP

#include "mymath_config.h"

#ifdef MYMATH_DBG
#include <string>
#include <sstream>
#endif

#include <Eigen/Dense>

#include "matrix_operations.hpp"

namespace mypilot {
namespace mymath {

/*
 * 离散实现的卡尔曼滤波器
 * 
 * XN 状态的维度
 * ZN 观察值的维度
 * UN 控制维度
 */
template <typename T, unsigned int XN, unsigned int ZN, unsigned int UN>
class KalmanFilter {
public:

  // 推迟初始化直到构造初始状态的构造函数,设置分布参数（使用SetStateEstimate）。
  // 通常在第一次观察时。
  KalmanFilter() {
    _F.setIdentity();
    _Q.setZero();
    _H.setIdentity();
    _R.setZero();
    _B.setZero();

    _x.setZero();
    _P.setZero();
    _y.setZero();
    _S.setZero();
    _K.setZero();
  }

  /*
   * x : 信念分布的期望。
   * P : 状态信念分布的协方差。
   */
  KalmanFilter(const Eigen::Matrix<T, XN, 1> &x,
               const Eigen::Matrix<T, XN, XN> &P)
      : KalmanFilter() {
    set_state_estimate(x, P);
  }

  virtual ~KalmanFilter() {}

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
    _is_initialized = true;
  }

  // 在零控制下更改系统转换功能
  void set_transition_matrix(const Eigen::Matrix<T, XN, XN>& F) { _F = F; }
  // 更改过渡噪声的协方差矩阵。
  void set_transition_noise(const Eigen::Matrix<T, XN, XN>& Q) { _Q = Q; }
  // 更改观察矩阵，将状态映射到观察。
  void set_observation_matrix(const Eigen::Matrix<T, ZN, XN>& H) { _H = H; }
  // 更改观测噪声的协方差矩阵。
  void set_observation_noise(const Eigen::Matrix<T, ZN, ZN>& R) { _R = R; }
  // 更改当前状态信念分布的协方差矩阵。
  void set_state_covariance(const Eigen::Matrix<T, XN, XN>& P) { _P = P; }
  // 更改状态转换规则中的控制矩阵。
  void set_control_matrix(const Eigen::Matrix<T, XN, UN>& B) { _B = B; }

  // 在零控制下获得系统转换函数。
  const Eigen::Matrix<T, XN, XN>& get_transition_matrix() const { return _F; }
  // 获取变换噪声的协方差矩阵。
  const Eigen::Matrix<T, XN, XN>& get_transition_noise() const { return _Q; }
  // 获取观测噪声的协方差矩阵。
  const Eigen::Matrix<T, ZN, XN>& get_observation_matrix() const { return _H; }
  // 获取观测噪声的协方差矩阵。
  const Eigen::Matrix<T, ZN, ZN>& get_observation_noise() const { return _R; }
  // 获取状态转换规则中的控制矩阵。
  const Eigen::Matrix<T, XN, UN>& get_control_matrix() const { return _B; }
  // 获取我们当前状态信念分布的平均值。
  Eigen::Matrix<T, XN, 1> get_state_estimate() const { return _x; }
  // 获取我们当前状态信念分布的协方差
  Eigen::Matrix<T, XN, XN> get_state_covariance() const { return _P; }

  // 是否初始化
  bool is_initialized() const { return _is_initialized; }

  // 给定控制输入'u'(默认为0)更新状态置信度分布。
  void predict(const Eigen::Matrix<T, UN, 1>& u=Eigen::Matrix<T, UN, 1>::Zero()) {
    assert(_is_initialized==true);
    _x = _F * _x + _B * u;
    _P = _F * _P * _F.transpose() + _Q;
  }

  // 给定观察值'z'更新状态置信分布。
  void correct(const Eigen::Matrix<T, ZN, 1>& z) {
    assert(_is_initialized==true);
    _y = z - _H * _x;
    _S = _H * _P * _H.transpose() + _R;
    _K = _P * _H.transpose() * pseudo_inverse<T, ZN>(_S);
    _x = _x + _K * _y;
    _P = (Eigen::Matrix<T, XN, XN>::Identity() - _K * _H) * _P;
  }

#if MYMATH_DBG
  std::string str() const {
    Eigen::IOFormat clean_fmt(4, 0, ", ", " ", "[", "]");
    std::ostringstream ss;
    ss << "F = " << _F.format(clean_fmt) << "\n"
      << "B = " << _B.format(clean_fmt) << "\n"
      << "H = " << _H.format(clean_fmt) << "\n"
      << "Q = " << _Q.format(clean_fmt) << "\n"
      << "R = " << _R.format(clean_fmt) << "\n"
      << "x = " << _x.format(clean_fmt) << "\n"
      << "P = " << _P.format(clean_fmt) << "\n";
    return ss.str();
  }
#endif

private:
  // 当前状态信念分布的期望
  Eigen::Matrix<T, XN, 1> _x;
  // 当前状态信念分布的协方差
  Eigen::Matrix<T, XN, XN> _P;
  // 零控制下的状态转移矩阵
  Eigen::Matrix<T, XN, XN> _F;
  // 状态转换噪声的协方差
  Eigen::Matrix<T, XN, XN> _Q;
  // 观察值矩阵
  Eigen::Matrix<T, ZN, XN> _H;
  // 观测噪声的协方差
  Eigen::Matrix<T, ZN, ZN> _R;
  // 状态转移规则中的控制矩阵
  Eigen::Matrix<T, XN, UN> _B;

  // 更新;标记为成员，以防止重新分配内存。
  Eigen::Matrix<T, ZN, 1> _y;
  // 更新协方差；标记为成员，以防止重新分配内存。
  Eigen::Matrix<T, ZN, ZN> _S;
  // 卡尔曼增益；标记为成员，以防止重新分配内存。
  Eigen::Matrix<T, XN, ZN> _K;

  // set_state_estimate被调用后设置为true。
  bool _is_initialized = false;
};

}}

#endif
