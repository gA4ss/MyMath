#ifndef MYMATH_LINEAR_QUADRATIC_REGULATOR_HPP
#define MYMATH_LINEAR_QUADRATIC_REGULATOR_HPP

#include <limits>

#include <Eigen/Core>
#include <Eigen/Dense>

namespace mypilot {
namespace mymath {

/**
 * @brief Solver for discrete-time linear quadratic problem.
 * @param A The system dynamic matrix
 * @param B The control matrix
 * @param Q The cost matrix for system state
 * @param R The cost matrix for control output
 * @param tolerance The numerical tolerance for solving
 *        Algebraic Riccati equation (ARE)
 * @param max_num_iteration The maximum iterations for solving ARE
 * @param ptr_K The feedback control matrix (pointer)
 */
void solve_lqr_problem(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B,
                       const Eigen::MatrixXd &Q, const Eigen::MatrixXd &R,
                       const double tolerance, const uint max_num_iteration,
                       Eigen::MatrixXd *ptr_K) {
  if (A.rows() != A.cols() || B.rows() != A.rows() || Q.rows() != Q.cols() ||
      Q.rows() != A.rows() || R.rows() != R.cols() || R.rows() != B.cols()) {
    // LQR solver: one or more matrices have incompatible dimensions.
    return;
  }

  Matrix AT = A.transpose();
  Matrix BT = B.transpose();

  // Solves a discrete-time Algebraic Riccati equation (DARE)
  // Calculate Matrix Difference Riccati Equation, initialize P and Q
  Matrix P = Q;
  uint num_iteration = 0;
  double diff = std::numeric_limits<double>::max();
  while (num_iteration++ < max_num_iteration && diff > tolerance) {
    Matrix P_next =
      AT * P * A - AT * P * B * (R + BT * P * B).inverse() * BT * P * A + Q;
    // check the difference between P and P_next
    diff = fabs((P_next - P).maxCoeff());
    P = P_next;
  }

  if (num_iteration >= max_num_iteration) {
    // LQR solver cannot converge to a solution, last consecutive result diff. is: diff
  } else {
    // LQR solver converged at iteration: num_iteration.
    // max consecutive result diff.: diff.
  }
  *ptr_K = (R + BT * P * B).inverse() * BT * P * A;
}

}}

#endif
