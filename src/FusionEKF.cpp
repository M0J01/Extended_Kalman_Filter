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


	H_laser_ << 1, 0, 0, 0,
							0, 1, 0, 0;

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

	// Maybe define P, Q, R?
	// ekf_.P_ << 1;
	// ekf_.Q_ << 1;
	// ekf_.R_ << 1;

	ekf_.P_ = MatrixXd(4,4);
	ekf_.P_ << 	1, 0, 1, 0,
							0, 1, 0, 1,
							0, 0, 1, 0,
							0, 0, 0, 1;

	//the initial transition matrix F_
	ekf_.F_ = MatrixXd(4, 4);
	ekf_.F_ << 1, 0, 1, 0,
					0, 1, 0, 1,
					0, 0, 1, 0,
					0, 0, 0, 1;

	ekf_.Q_ = MatrixXd(4, 4);
	ekf_.Q_ << 	1, 0, 1, 0,
							0, 1, 0, 1,
							1, 0, 1, 0,
							0, 1, 0, 1;
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
		previous_timestamp_ = measurement_pack.timestamp_;

    // First Measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

		// ---- RADAR ----
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

			// Convert Polar to Cartesian
			float ro = measurement_pack.raw_measurements_(0);
			float phi = measurement_pack.raw_measurements_(1);
			//float ro_dot = measurement_pack.raw_measurements_(2); // Not Enough Info for this to work
			float p_x = ro * cos(phi);
			float p_y = ro * sin(phi);

			// Load position data in our state matrix x_
			ekf_.x_ << p_x, p_y, 1, 1;

		}
    // ---- LIDAR ----
		else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

			// Load data in our state matrix x_
			float x = measurement_pack.raw_measurements_(0);
			float y = measurement_pack.raw_measurements_(1);
			ekf_.x_ << x, y, 1, 1;

    }

    // done initializing, no need to predict or update

    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/


	long double dt = (measurement_pack.timestamp_ - previous_timestamp_ ) /1000000.0 ;
	previous_timestamp_ = measurement_pack.timestamp_;


	// Update our state transition matrix F_ based on our time elapsed
	ekf_.F_ << 	1, 0, dt, 0,
							0, 1, 0, dt,
							0, 0, 1, 0,
							0, 0, 0, 1;

	// Set Q_ parameters
	float noise_ax = 81;
	float noise_ay = 81;
	double ax4 = pow(dt,4)*noise_ax/4.0;
	double ax3 = pow(dt,3)*noise_ax/2.0;
	double ax2 = pow(dt,2)*noise_ax;
	double ay4 = pow(dt,4)*noise_ay/4.0;
	double ay3 = pow(dt,3)*noise_ay/2.0;
	double ay2 = pow(dt,2)*noise_ay;

	// Update our Covariance Q_ based on our time elapsed

	ekf_.Q_ << 	ax4, 0, ax3, 0,
							0, ay4, 0, ay3,
							ax3, 0, ax2, 0,
							0, ay3, 0, ay2;

	// Predict new state x and process covariance P matrices.
	ekf_.Predict();


  /*****************************************************************************
   *  Update
   ****************************************************************************/


	// Radar updates
	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {


		// Calculate our Hj_ matrix

		// cout << "This is what Hj looks like.\n-->\n";
		Hj_ = tools.CalculateJacobian(ekf_.x_);
		//Hj_ = tools.CalculateJacobian(measurement_pack.raw_measurements_);
		 //cout << Hj_ << "\n<--\nThis is what Hj looks like.";

		//cout << "Radar step 2!";

		// Set our H and R matrices
		ekf_.H_ = Hj_;
		ekf_.R_ = R_radar_;

		// Perform Update step
		ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {

		// Set our H and R matrices
		ekf_.H_ = H_laser_;
		ekf_.R_ = R_laser_;

		// Perform our Lidar Update step
		ekf_.Update(measurement_pack.raw_measurements_);
	}

  // print the output
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;

 }
