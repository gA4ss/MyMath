#ifndef MYMATH_QP_SOLVER_ACTIVE_SET_QP_SOLVER_HPP
#define MYMATH_QP_SOLVER_ACTIVE_SET_QP_SOLVER_HPP

#include <cassert>
#include <climits>
#include <vector>
#include <algorithm>
#include <qpOASES.hpp>

#include "qp_solver_gflags.h"
#include "qp_solver.hpp"

namespace mypilot {
namespace mymath {

class ActiveSetQpSolver : public QpSolver {
public:
  ActiveSetQpSolver(const Eigen::MatrixXd& kernel_matrix,
                    const Eigen::MatrixXd& offset,
                    const Eigen::MatrixXd& affine_inequality_matrix,
                    const Eigen::MatrixXd& affine_inequality_boundary,
                    const Eigen::MatrixXd& affine_equality_matrix,
                    const Eigen::MatrixXd& affine_equality_boundary) : 
    QpSolver(kernel_matrix, offset, affine_inequality_matrix,
             affine_inequality_boundary, affine_equality_matrix,
             affine_equality_boundary),
    _num_constraint(_affine_equality_matrix.rows() +
                    _affine_inequality_matrix.rows()),
    _num_param(kernel_matrix.rows()),
    _qp_eps_num(default_active_set_eps_num),
    _qp_eps_den(default_active_set_eps_den),
    _qp_eps_iter_ref(default_active_set_eps_iter_ref),
    _debug_info(default_enable_active_set_debug_info) {}
  virtual ~ActiveSetQpSolver() = default;

  bool solve() override {
    ::qpOASES::QProblem qp_problem(_num_param, _num_constraint, _hessian_type);
    ::qpOASES::Options my_options;

    my_options.enableCholeskyRefactorisation = _cholesky_refactorisation_freq;
    if (_hessian_type == ::qpOASES::HST_POSDEF ||
        _hessian_type == ::qpOASES::HST_SEMIDEF) {
      my_options.enableRegularisation = ::qpOASES::BT_TRUE;
    }
    my_options.epsNum = _qp_eps_num;
    my_options.epsDen = _qp_eps_den;
    my_options.epsIterRef = _qp_eps_iter_ref;
    my_options.terminationTolerance = _termination_tolerance;
    qp_problem.setOptions(my_options);
    if (!_debug_info) {
      qp_problem.setPrintLevel(qpOASES::PL_NONE);
    }
    if (_kernel_matrix.rows() != _kernel_matrix.cols()) {
      /*
        stringstream ss; 
        ss << "_kernel_matrix.rows() [" << _kernel_matrix.rows()
          << "] and _kernel_matrix.cols() [" << _kernel_matrix.cols()
          << "] should be identical.";
      */
      return false;
    }
    // 定义qpOASESproblem
    const int kNumOfMatrixElements =
      _kernel_matrix.rows() * _kernel_matrix.cols();
    double h_matrix[kNumOfMatrixElements];

    const int kNumOfOffsetRows = _offset.rows();
    double g_matrix[kNumOfOffsetRows];
    int index = 0;

    for (int r = 0; r < _kernel_matrix.rows(); ++r) {
      g_matrix[r] = _offset(r, 0);
      for (int c = 0; c < _kernel_matrix.cols(); ++c) {
        h_matrix[index++] = _kernel_matrix(r, c);
      }
    }
    assert(index == _kernel_matrix.rows() * _kernel_matrix.cols());

    // 搜索空间的下限与上限
    double lower_bound[_num_param];
    double upper_bound[_num_param];

    // 配置上下限
    for (int i = 0; i < _num_param; ++i) {
      lower_bound[i] = _l_lower_bound;
      upper_bound[i] = _l_upper_bound;
    }

    // 约束矩阵构造
    double affine_constraint_matrix[_num_param * _num_constraint];
    double constraint_lower_bound[_num_constraint];
    double constraint_upper_bound[_num_constraint];
    index = 0;

    for (int r = 0; r < _affine_equality_matrix.rows(); ++r) {
      // 配置上下限
      constraint_lower_bound[r] = _affine_equality_boundary(r, 0);
      constraint_upper_bound[r] = _affine_equality_boundary(r, 0);

      for (int c = 0; c < _num_param; ++c) {
        affine_constraint_matrix[index++] = _affine_equality_matrix(r, c);
      }
    }

    assert(index == _affine_equality_matrix.rows() * _num_param);

    for (int r = 0; r < _affine_inequality_matrix.rows(); ++r) {
      constraint_lower_bound[r + _affine_equality_boundary.rows()] =
          _affine_inequality_boundary(r, 0);
      constraint_upper_bound[r + _affine_equality_boundary.rows()] =
          _constraint_upper_bound;

      // 配置上下限
      for (int c = 0; c < _num_param; ++c) {
        affine_constraint_matrix[index++] = _affine_inequality_matrix(r, c);
      }
    }
    assert(index == _affine_equality_matrix.rows() * _num_param +
                    _affine_inequality_boundary.rows() * _num_param);

    // 初始化
    int max_iter = std::max(_max_iteration, _num_constraint);

    auto ret = qp_problem.init(h_matrix, g_matrix, affine_constraint_matrix,
                               lower_bound, upper_bound, constraint_lower_bound,
                               constraint_upper_bound, max_iter);

    //
    // 解决问题失败
    //
    if (ret != qpOASES::SUCCESSFUL_RETURN) {
      if (ret == qpOASES::RET_MAX_NWSR_REACHED) {
        // qpOASES solver failed due to reached max iteration
      } else {
        // qpOASES solver failed due to infeasibility or other internal reasons
      }
      std::stringstream ss;
      ss << "ActiveSetQpSolver inputs: " << std::endl;
      ss << "kernel_matrix:\n" << _kernel_matrix << std::endl;
      ss << "offset:\n" << _offset << std::endl;
      ss << "affine_inequality_matrix:\n"
         << _affine_inequality_matrix << std::endl;
      ss << "affine_inequality_boundary:\n"
         << _affine_inequality_boundary << std::endl;
      ss << "affine_equality_matrix:\n" << _affine_equality_matrix << std::endl;
      ss << "affine_equality_boundary:\n"
         << _affine_equality_boundary << std::endl;
      // ss.str();
      return false;
    }

    double result[_num_param];  // NOLINT
    qp_problem.getPrimalSolution(result);

    _params = Eigen::MatrixXd::Zero(_num_param, 1);
    for (int i = 0; i < _num_param; ++i) {
      _params(i, 0) = result[i];
    }
    return qp_problem.isSolved() == qpOASES::BT_TRUE;
  }

  void set_qp_eps_num(const double eps) { _qp_eps_num = eps; }
  void set_qp_eps_den(const double eps) { _qp_eps_den = eps; }
  void set_qp_eps_iter_ref(const double eps) { _qp_eps_iter_ref = eps; }
  void set_debug_info(const bool enable) { _debug_info = enable; }
  void set_max_iteration(const int max_iter) { _max_iteration = max_iter; }

  void set_l_lower_bound(const double l_lower_bound) { _l_lower_bound = l_lower_bound; }
  void set_l_upper_bound(const double l_upper_bound) { _l_upper_bound = l_upper_bound; }
  void set_constraint_upper_bound(const double la_upper_bound) { _constraint_upper_bound = la_upper_bound; }

  double qp_eps_num() const { return _qp_eps_num; }
  double qp_eps_den() const { return _qp_eps_den; }
  double qp_eps_iter_ref() const { return _qp_eps_iter_ref; }
  bool debug_info() const { return _debug_info; }
  int max_iteration() const { return _max_iteration; }

  double l_lower_bound() const { return _l_lower_bound; }
  double l_upper_bound() const { return _l_upper_bound; }
  double constraint_upper_bound() const { return _constraint_upper_bound; }

  void set_pos_semi_definite_hessian() override {
    // 如果'hessian'类型设置为'HST_SEMIDEF'
    // 内置的正则化方案无需额外的计算就可以打开
    _hessian_type = ::qpOASES::HST_SEMIDEF;
  }
  void set_pos_definite_hessian() override {
    _hessian_type = ::qpOASES::HST_POSDEF;
  }

  void enable_cholesky_refactorisation(const int num) override {
    // 指定预计的Hessian的完全分解重构的频率：0将其关闭，1在每次迭代中使用它们，等等。
    _cholesky_refactorisation_freq = num;
  }

  void set_termination_to_lerance(const double tolerance) override {
    _termination_tolerance = tolerance;
  }

private:
  bool sanity_check() override {
    return _kernel_matrix.rows() == _kernel_matrix.cols() &&
          _kernel_matrix.rows() == _affine_inequality_matrix.cols() &&
          _kernel_matrix.rows() == _affine_equality_matrix.cols() &&
          _affine_equality_matrix.rows() == _affine_equality_boundary.rows() &&
          _affine_inequality_matrix.rows() == _affine_inequality_boundary.rows();
  }

private:
  int _num_constraint = 0;      // 等式约束数量 + 不等式约束数量
  int _num_param = 0;           // 参数数量

  double _qp_eps_num = 0.0;
  double _qp_eps_den = 0.0;
  double _qp_eps_iter_ref = 0.0;
  bool _debug_info = false;

  // 搜索界限参数
  double _l_lower_bound = -1e10;
  double _l_upper_bound = 1e10;

  // 搜索上界约束
  double _constraint_upper_bound = 1e10;
  int _max_iteration = 1000;

  ::qpOASES::HessianType _hessian_type = ::qpOASES::HST_UNKNOWN;
  int _cholesky_refactorisation_freq = 0;
  double _termination_tolerance = 1.0e-9;
};

}}

#endif
