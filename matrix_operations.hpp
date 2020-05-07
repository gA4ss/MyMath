#ifndef MYMATH_MATRIX_OPERATIONS_HPP
#define MYMATH_MATRIX_OPERATIONS_HPP

#include <cmath>
#include <utility>
#include <cassert>

#include <Eigen/Dense>
#include <Eigen/SVD>

namespace mypilot {
namespace mymath {

// 计算给定方矩阵的Moore-Penrose伪逆，将具有以epsilon为界的绝对值的所有特征值四舍五入。
template <typename T, unsigned int N>
Eigen::Matrix<T, N, N> pseudo_inverse(const Eigen::Matrix<T, N, N>& m,
                                      const double epsilon=1.0e-6) {
  // Jacobi的特征值求解算法
  Eigen::JacobiSVD<Eigen::Matrix<T, N, N>> svd(
    m, Eigen::ComputeFullU | Eigen::ComputeFullV);
  return svd.matrixV() *
    (svd.singularValues().array().abs() > epsilon)
        .select(svd.singularValues().array().inverse(), 0)
        .matrix()
        .asDiagonal() *
    svd.matrixU().adjoint();
}

// 计算给定矩阵的Moore-Penrose伪逆，将具有以epsilon为界的绝对值的所有特征值四舍五入。
template <typename T, unsigned int M, unsigned int N>
Eigen::Matrix<T, N, M> pseudo_inverse(const Eigen::Matrix<T, M, N>& m,
                                      const double epsilon=1.0e-6) {
  Eigen::Matrix<T, M, M> t = m * m.transpose();
  return m.transpose() * pseudo_inverse<T, M>(t);
}

/*
 * 计算状态空间表示形式的连续形式到离散形式的双线性变换。假设方程式为 ：
 * 
 *         dot_x = Ax + Bu
 *         y = Cx + Du
 * 参数 m_a, m_b, m_c, m_d 是状态空间的控制矩阵。
 */
template <typename T, unsigned int L, unsigned int N, unsigned int O>
bool continuous_to_discrete(const Eigen::Matrix<T, L, L>& m_a,
                            const Eigen::Matrix<T, L, N>& m_b,
                            const Eigen::Matrix<T, O, L>& m_c,
                            const Eigen::Matrix<T, O, N>& m_d, const double ts,
                            Eigen::Matrix<T, L, L>* ptr_a_d,
                            Eigen::Matrix<T, L, N>* ptr_b_d,
                            Eigen::Matrix<T, O, L>* ptr_c_d,
                            Eigen::Matrix<T, O, N>* ptr_d_d) {
  if (ts <= 0.0) {
    // continuous_to_discrete : ts is less than or equal to zero
    return false;
  }

  // 在矩阵中，只有matrix_a强制为非零。
  if (m_a.rows() == 0) {
    // continuous_to_discrete: matrix_a size 0
    return false;
  }

  Eigen::Matrix<T, L, L> m_identity = Eigen::Matrix<T, L, L>::Identity();
  *ptr_a_d = pseudo_inverse<T, L>(m_identity - ts * 0.5 * m_a) *
    (m_identity + ts * 0.5 * m_a);
  *ptr_b_d =
    std::sqrt(ts) * pseudo_inverse<T, L>(m_identity - ts * 0.5 * m_a) * m_b;
  *ptr_c_d =
    std::sqrt(ts) * m_c * pseudo_inverse<T, L>(m_identity - ts * 0.5 * m_a);
  *ptr_d_d =
    0.5 * m_c * pseudo_inverse<T, L>(m_identity - ts * 0.5 * m_a) * m_b + m_d;
  return true;
}

bool continuous_to_discrete(const Eigen::MatrixXd& m_a,
                            const Eigen::MatrixXd& m_b,
                            const Eigen::MatrixXd& m_c,
                            const Eigen::MatrixXd& m_d, const double ts,
                            Eigen::MatrixXd* ptr_a_d, Eigen::MatrixXd* ptr_b_d,
                            Eigen::MatrixXd* ptr_c_d, Eigen::MatrixXd* ptr_d_d) {
  if (ts <= 0.0) {
    // continuous_to_discrete : ts is less than or equal to zero
    return false;
  }

  // 在矩阵中，只有matrix_a强制为非零。
  if (m_a.rows() == 0) {
    // continuous_to_discrete : matrix_a size 0
    return false;
  }

  if (m_a.cols() != m_b.rows() || m_b.cols() != m_d.cols() ||
      m_c.rows() != m_d.rows() || m_a.cols() != m_c.cols()) {
    // continuous_to_discrete: matrix dimensions mismatch
    return false;
  }

  Eigen::MatrixXd m_identity =
    Eigen::MatrixXd::Identity(m_a.cols(), m_a.rows());

  *ptr_a_d =
    (m_identity - ts * 0.5 * m_a).inverse() * (m_identity + ts * 0.5 * m_a);
  *ptr_b_d = std::sqrt(ts) * (m_identity - ts * 0.5 * m_a).inverse() * m_b;
  *ptr_c_d = std::sqrt(ts) * m_c * (m_identity - ts * 0.5 * m_a).inverse();
  *ptr_d_d = 0.5 * m_c * (m_identity - ts * 0.5 * m_a).inverse() * m_b + m_d;
  return true;
}

}}

#endif
