#ifndef MYMATH_QP_SOLVER_QP_SOLVER_GFLAGS_H
#define MYMATH_QP_SOLVER_QP_SOLVER_GFLAGS_H

namespace mypilot {
namespace mymath {

// qpOases wrapper error control numerator
#define default_active_set_eps_num              -1e-7
// qpOases wrapper error control denominator
#define default_active_set_eps_den              1e-7
// qpOases wrapper early termination tolerance for iterative refinement
#define default_active_set_eps_iter_ref         1e-7
// qpOases wrapper error control numerator
#define default_qp_smoothing_eps_num            -1e-7
// qpOases wrapper error control denominator
#define default_qp_smoothing_eps_den            1e-7
// qpOases wrapper early termination tolerance for iterative refinement
#define default_qp_smoothing_eps_iter_ref       1e-7
// Enable print information
#define default_enable_active_set_debug_info    false
// Default qp oases iteration time
#define default_qp_iteration_num                1000

}}

#endif