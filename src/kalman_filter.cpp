#include "kalman_filter.h"
#include <math.h>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
		x_ = F_ * x_;
		P_ = F_ * P_ * F_.transpose() + Q_;

		/*std::cout << "PREDICT" << std::endl;
		std::cout << "x=" << std::endl << x_ << std::endl;
		std::cout << "P=" << std::endl << P_ << std::endl;
		std::cout << "END PREDICT" << std::endl;*/
}

void KalmanFilter::Update(const VectorXd &z) {  
		VectorXd y = z - (H_ * x_);
		MatrixXd S = H_ * P_ * H_.transpose() + R_;
		MatrixXd K = P_ * H_.transpose() * S.inverse();
		x_ = x_ + (K * y);
		long x_size = x_.size();
		MatrixXd I = MatrixXd::Identity(x_size, x_size);
		P_ = (I - (K * H_)) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {	
		VectorXd h_x(3);	
		double c1 = sqrt(x_[0]*x_[0] + x_[1]*x_[1]);
		double c2 = (x_[1] != 0 && x_[0] != 0) ? atan2(x_[1], x_[0]) : 0;
		const double epsilon = 0.00001;
		h_x << c1, c2, (x_[0]*x_[2] + x_[1]*x_[3])/std::max(c1, epsilon);
		VectorXd y = z - h_x;
		KalmanFilter::NormaliseAngle(y(1));
		MatrixXd Hj = tools_.CalculateJacobian(x_);
		MatrixXd S = Hj * P_ * Hj.transpose() + R_;
		MatrixXd K = P_ * Hj.transpose() * S.inverse();
		x_ = x_ + (K * y);
		long x_size = x_.size();
		MatrixXd I = MatrixXd::Identity(x_size, x_size);
		P_ = (I - (K * Hj)) * P_;
		/*std::cout << "UPDATE" << std::endl;
		std::cout << "y=" << std::endl << y << std::endl;
		std::cout << "Hj=" << std::endl << Hj << std::endl;
		std::cout << "S=" << std::endl << S << std::endl;
		std::cout << "K=" << std::endl << K << std::endl;
		std::cout << "x=" << std::endl << x_ << std::endl;
		std::cout << "P=" << std::endl << P_ << std::endl;
		std::cout << "END UPDATE" << std::endl;*/
}

void KalmanFilter::NormaliseAngle(double &angle) {
	while (angle < -M_PI || angle > M_PI) {
		double increment = angle > 0 ? -2 * M_PI : 2 * M_PI;
		angle += increment;
	}
}
