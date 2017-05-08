#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd y = z - H_ * x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();
  x_ = x_ + K * y;
  MatrixXd I(4, 4);
  I << 1., 0, 0, 0,
    0, 1., 0, 0,
    0, 0, 1., 0,
    0, 0, 0, 1.;
  P_ = (I - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  VectorXd hx(3);
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  hx(0) = sqrt(px * px + py * py);
  hx(1) = atan2(py, px);
  hx(2) = (vy * py + vx * px)/hx(0);
  VectorXd y = z - hx;
  if (y(1) > M_PI) {
    y(1) = y(1) - 2 * M_PI;
  }
  else if (y(1) < -M_PI) {
    y(1) = y(1) + 2 * M_PI;
  }
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();
  x_ = x_ + K * y;
  MatrixXd I(4, 4);
  I << 1., 0, 0, 0,
    0, 1., 0, 0,
    0, 0, 1., 0,
    0, 0, 0, 1.;
  P_ = (I - K * H_) * P_;
}
