#include "kalman_filter.h"
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

	// Perform Prediction
	x_ = F_ * x_; // + u_, which we set to 0, as it shows up in our Q Process Covariance
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;

	// std::cout << "\n\n This is X during predict: x = \n" << x_;
}

void KalmanFilter::Update(const VectorXd &z) {

	long x_size = x_.size();
	MatrixXd I_ = MatrixXd::Identity(x_size, x_size);

	// Perform Update
	VectorXd y = z - H_* x_;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K = P_ * Ht * Si;
	x_ = x_ + (K * y);
	P_ = (I_ - K * H_) * P_;

	// std::cout << "\n\n This is X during Update: x = \n" << x_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

	long x_size = x_.size();
	MatrixXd I_ = MatrixXd::Identity(x_size, x_size);

	// Compute h(x) for our radar
	double p_x = x_(0);
	double p_y = x_(1);
	double v_x  = x_(2);
	double v_y = x_(3);
	//std::cout << std::endl << " Px = " << p_x << ", Py = " << p_y << ", Vx = " << v_x << ", Vy = " << v_y << std::endl;
	// Catch divide by 0
	float denom = pow(p_x, 2) + pow(p_y, 2);
	if (fabs(denom) < 0.00000001) { denom = 0.00000001; }

	// Convert Cartesian to Polar
	float ro = sqrt(pow(p_x, 2) + pow(p_y, 2));

	if (fabs(p_x) < 0.0001) { p_x = 0.0001;}
	float phi = (atan2(p_y,p_x));

	//if (phi < 0){phi = -3.14 - phi;}

	float ro_dot = ((p_x*v_x + p_y*v_y)/sqrt(denom));
	if (z(1) > 0 && phi < 0){ phi = fabs(phi);}

	// std::cout << "\n\nro = " << ro << ", phi = " << phi << ", ro' = " << ro_dot;

	// Load our h(x) vector
	VectorXd h_x(3); // (p, phi, p*)
	h_x << 	ro,
					phi,
					ro_dot;

	// Perform Update
	VectorXd y = z - h_x;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K = P_ * Ht * Si;
	x_ = x_ + (K * y);
	P_ = (I_ - K * H_) * P_;

	//std::cout << "\n z = \n" << z;
	//std::cout << "\ny = \n" << y;
	//std::cout << "\n H = \n" << H_;
	//std::cout<< "Ht = \n" << Ht;
	//std::cout << "\n S = \n" << S;
	//std::cout << "\n Si = \n" << Si;
	//std::cout << "\n K = \n" << K;
	//std::cout << "\n P = \n" << P_;

	// std::cout << "\n\n This is X during Update: x = \n" << x_;


	//std::getchar();
}
