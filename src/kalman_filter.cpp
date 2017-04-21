#include "kalman_filter.h"
#include <math.h>

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
		float c1 = sqrt(x_[0]*x_[0] + x_[1]*x_[1]);
		float c2 = atan(_x[1]/x_[0]);
		while(c2 < -PI || c2 > PI){
			float increment = c2 > 0 ? -2*PI : 2*PI;
			c2 += increment;
		}
		h_x << c1, c2, (x_[0]*x_[2] + x_[1]*x_[4])/c1;
		VectorXd y = z - h_x;
		Hj = tools_.CalculateJacobian(x_);
		MatrixXd S = Hj * P_ * Hj.transpose() + R_;
		MatrixXd K = P_ * Hj.transpose() * S.inverse();
		x_ = x_ + (K * y);
		long x_size = x_.size();
		MatrixXd I = MatrixXd::Identity(x_size, x_size);
		P_ = (I - (K * Hj)) * P_;
}
