#ifndef MYMATH_QP_SOLVER_QP_SOLVER_HPP
#define MYMATH_QP_SOLVER_QP_SOLVER_HPP

#include <Eigen/Core>
#include <Eigen/LU>

namespace mypilot {
namespace mymath {

/*
 * 二次多项式问题求解
 * 
 * min_x  : q(x) = 0.5 * x^T * Q * x  + x^T c
 * A * x = b (等式约束)
 * C * x >= d (不等式约束)
 */
class QpSolver {
public:
  QpSolver(const Eigen::MatrixXd& kernel_matrix, const Eigen::MatrixXd& offset,
           const Eigen::MatrixXd& affine_inequality_matrix,
           const Eigen::MatrixXd& affine_inequality_boundary,
           const Eigen::MatrixXd& affine_equality_matrix,
           const Eigen::MatrixXd& affine_equality_boundary) : 
    _kernel_matrix(kernel_matrix),
    _offset(offset),
    _affine_inequality_matrix(affine_inequality_matrix),
    _affine_inequality_boundary(affine_inequality_boundary),
    _affine_equality_matrix(affine_equality_matrix),
    _affine_equality_boundary(affine_equality_boundary) {}

  virtual ~QpSolver() = default;

  virtual void set_pos_semi_definite_hessian() {}
  virtual void set_pos_definite_hessian() {}
  virtual void enable_cholesky_refactorisation(const int) {}
  virtual void set_termination_to_lerance(const double) {}
  virtual bool solve() = 0;

  const Eigen::MatrixXd& params() const { return _params; }
  const Eigen::MatrixXd& kernel_matrix() const { return _kernel_matrix; }
  const Eigen::MatrixXd& offset() const { return _offset; }
  const Eigen::MatrixXd& affine_equality_matrix() const { return _affine_equality_matrix; }
  const Eigen::MatrixXd& affine_equality_boundary() const { return _affine_equality_boundary; }
  const Eigen::MatrixXd& affine_inequality_matrix() const { return _affine_inequality_matrix; }
  const Eigen::MatrixXd& affine_inequality_boundary() const { return _affine_inequality_boundary; }

protected:
  virtual bool sanity_check() = 0;
  Eigen::MatrixXd _params;
  Eigen::MatrixXd _kernel_matrix;
  Eigen::MatrixXd _offset;
  Eigen::MatrixXd _affine_inequality_matrix;
  Eigen::MatrixXd _affine_inequality_boundary;
  Eigen::MatrixXd _affine_equality_matrix;
  Eigen::MatrixXd _affine_equality_boundary;
};

}}

#endif
