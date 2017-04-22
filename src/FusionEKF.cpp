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

  H_laser_ << 1, 0, 0, 0,
	  0, 1, 0, 0;

  Hj_ << 0, 0, 0, 0,
	  0, 0, 0, 0,
	  0, 0, 0, 0;


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
    // first measurement
    //cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
		double rho = measurement_pack.raw_measurements_[0];
		double phi = measurement_pack.raw_measurements_[1];
		ekf_.x_ << rho * cos(phi), rho * sin(phi), 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
		ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

	//the initial transition matrix F_
	ekf_.F_ = MatrixXd(4, 4);
	ekf_.F_ << 1, 0, 1, 0,
		0, 1, 0, 1,
		0, 0, 1, 0,
		0, 0, 0, 1;

	//state covariance matrix P
	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1000, 0,
		0, 0, 0, 1000;


	//set the process covariance matrix Q
	ekf_.Q_ = MatrixXd(4, 4);

	previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  int noise_ax = 9;
  int noise_ay = 9;
  double dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  double pow_2 = dt * dt;
  double pow_3 = pow_2 * dt;
  double pow_4 = pow_3 * dt;
  double Q_2_4 = pow_3 * noise_ay / 2.0;
  double Q_1_3 = pow_3 * noise_ax / 2.0;
  ekf_.Q_ << pow_4 * noise_ax / 4.0, 0, Q_1_3, 0,
	  0, pow_4 * noise_ay / 4.0, 0, Q_2_4,
	  Q_1_3, 0, pow_2 * noise_ax, 0,
	  0, Q_2_4, 0, pow_2 * noise_ay;
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/


  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
	ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
	ekf_.H_ = H_laser_;
	ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}
