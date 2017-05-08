#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  // init laser
  H_laser_ << 1, 0, 0, 0,
    0, 1, 0, 0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float ro = measurement_pack.raw_measurements_(0);
      float theta = measurement_pack.raw_measurements_(1);
      float ro_dot = measurement_pack.raw_measurements_(2);
      ekf_.x_(0) = ro * sin(theta);
      ekf_.x_(1) = ro * cos(theta);
      // dubious: the following could be inaccurate but still best guess
      ekf_.x_(2) = ro_dot * sin(theta);
      ekf_.x_(3) = ro_dot * cos(theta);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_(0) =  measurement_pack.raw_measurements_(0);
      ekf_.x_(1) =  measurement_pack.raw_measurements_(1);
      ekf_.x_(2) =  0;
      ekf_.x_(3) =  0;
    }
    // initialize uncertainty matrix to have large uncertainty
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1000, 0,
      0, 0, 0, 1000;

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  float delta_t = (measurement_pack.timestamp_ - previous_timestamp_)/1e6;
  previous_timestamp_ = measurement_pack.timestamp_;
  // initialize Q
  MatrixXd sigma(2, 2);
  sigma << 9, 0,
    0, 9;
  MatrixXd G(4, 2);
  G << 0.5 * delta_t * delta_t, 0,
    0, 0.5 * delta_t * delta_t,
    delta_t, 0,
    0, delta_t;
  ekf_.Q_ = G * sigma * G.transpose();
  MatrixXd Q(4, 4);
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, delta_t, 0,
    0, 1, 0, delta_t,
    0, 0, 1, 0,
    0, 0, 0, 1;

  ekf_.Predict();
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    float ro = measurement_pack.raw_measurements_(0);
    float theta = measurement_pack.raw_measurements_(1);
    float ro_dot = measurement_pack.raw_measurements_(2);
    VectorXd z(3);
    z << ro, theta, ro_dot;
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(z);
  } else {
    // Laser updates
    float px = measurement_pack.raw_measurements_(0);
    float py = measurement_pack.raw_measurements_(1);
    VectorXd z(2);
    z << px, py;
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(z);
  }


  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
