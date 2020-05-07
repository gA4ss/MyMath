#include "kalman_filter.hpp"
#include "ltest.hpp"

using namespace mypilot::mymath;

class KalmanFilterTest {
public:
  KalmanFilterTest() : kf_() {}

  virtual void setup() {
    // Initial state
    Eigen::Matrix<double, 2, 1> x;
    x(0, 0) = 0.0;
    x(1, 0) = 1.0;

    // Initial state belief covariance
    Eigen::Matrix<double, 2, 2> P;
    P.setZero();
    P(0, 0) = 0.1;
    P(1, 1) = 0.1;

    // Transition matrix
    Eigen::Matrix<double, 2, 2> F;
    F.setZero();
    F(0, 0) = 1.0;
    F(0, 1) = 1.0;
    F(1, 1) = 1.0;

    // Transition noise covariance
    Eigen::Matrix<double, 2, 2> Q;
    Q.setZero();
    Q(0, 0) = 0.01;
    Q(1, 1) = 0.01;

    // Observation matrix
    Eigen::Matrix<double, 1, 2> H;
    H.setIdentity();

    // Observation noise covariance
    Eigen::Matrix<double, 1, 1> R;
    R(0, 0) = 0.5 * 0.5;

    // Control matrix
    Eigen::Matrix<double, 2, 1> B;
    B[0] = 0.5 * 1.0 * 1.0;
    B[1] = 1.0;

    kf_.set_state_estimate(x, P);
    kf_.set_transition_matrix(F);
    kf_.set_transition_noise(Q);
    kf_.set_observation_matrix(H);
    kf_.set_observation_noise(R);
    kf_.set_control_matrix(B);
  }

  void testing() {
    TEST_START("synthetic tracking test");
    {
      kf_.predict();
      Eigen::Matrix<double, 2, 1> state = kf_.get_state_estimate();
      Eigen::Matrix<double, 2, 2> state_cov = kf_.get_state_covariance();
      EXPECT_DOUBLE_EQ(1.0, state(0, 0));
      EXPECT_DOUBLE_EQ(1.0, state(1, 0));
      EXPECT_NEAR(0.21, state_cov(0, 0), 0.001);
      EXPECT_NEAR(0.10, state_cov(0, 1), 0.001);
      EXPECT_NEAR(0.10, state_cov(1, 0), 0.001);
      EXPECT_NEAR(0.11, state_cov(1, 1), 0.001);

      Eigen::Matrix<double, 1, 1> z;
      z(0, 0) = 1.0;
      kf_.correct(z);
      state = kf_.get_state_estimate();
      state_cov = kf_.get_state_covariance();

      EXPECT_DOUBLE_EQ(1.0, state(0, 0));
      EXPECT_DOUBLE_EQ(1.0, state(1, 0));
      EXPECT_NEAR(0.11413, state_cov(0, 0), 0.001);
      EXPECT_NEAR(0.05348, state_cov(0, 1), 0.001);
      EXPECT_NEAR(0.05348, state_cov(1, 0), 0.001);
      EXPECT_NEAR(0.08826, state_cov(1, 1), 0.001);
    }
    TEST_END("synthetic tracking test");
  }

protected:
  KalmanFilter<double, 2, 1, 1> kf_;
};

int main(int argc, char* argv[]) {
  KalmanFilterTest kf_test = KalmanFilterTest();
  kf_test.setup();
  kf_test.testing();
  return 0;
}
