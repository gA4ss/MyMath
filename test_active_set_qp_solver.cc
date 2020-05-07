#include "qp_solver/active_set_qp_solver.hpp"
#include "ltest.hpp"

using namespace mypilot::mymath;
using Eigen::MatrixXd;

int main(int argc, char* argv[]) {

  TEST_START("simple_problem_01");
  {
    MatrixXd kernel_matrix = MatrixXd::Zero(1, 1);
    kernel_matrix(0, 0) = 1.0;
    MatrixXd offset = MatrixXd::Zero(1, 1);
    offset(0, 0) = -8.0;
    MatrixXd affine_inequality_matrix;
    MatrixXd affine_inequality_boundary;
    MatrixXd affine_equality_matrix;
    MatrixXd affine_equality_boundary;
    ActiveSetQpSolver solver(kernel_matrix, offset, affine_inequality_matrix,
                            affine_inequality_boundary, affine_equality_matrix,
                            affine_equality_boundary);
    solver.solve();
    EXPECT_NEAR(solver.params()(0, 0), 8.0, 1e-9);
  }
  TEST_END("simple_problem_01");
}
