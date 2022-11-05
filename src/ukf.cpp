#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  // P_ = MatrixXd(5, 5);
  P_ = MatrixXd::Identity(5, 5);
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 2.5;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

  // State dimension
  n_x_ = 5; // todo

  // Augmented state dimension
  n_aug_ = 7; // todo

  // Sigma point spreading parameter
  lambda_ = 3 - n_x_; // todo
  lambda_aug_ = 3 - n_aug_;// 
  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_ + 1);

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1); // to do
  Xsig_pred_.fill(0.0);
  // compute weights
  double weight_0 = lambda_aug_/(lambda_aug_+n_aug_);
  double weight = 0.5/(lambda_aug_ + n_aug_);
  weights_(0) = weight_0;

  for (int i=1; i<2*n_aug_+1; ++i) {  
    weights_(i) = weight;
  }
  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false; // todo

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   * 
   */
    if (is_initialized_ == false){
      if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
        x_ << meas_package.raw_measurements_[0], 
              meas_package.raw_measurements_[1], 
              0, 
              0,
              0;
        P_(0,0) = std_laspx_*std_laspx_;
        P_(1,1) = std_laspy_*std_laspy_;  
      } 
      else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
        x_ << meas_package.raw_measurements_[0]*cos(meas_package.raw_measurements_[1]), 
              meas_package.raw_measurements_[0]*sin(meas_package.raw_measurements_[1]), 
              0, 
              0,
              0;
        P_ = P_ * std_radr_*std_radr_;
      }
      is_initialized_ = true;
    }
    else
    {
      double dt = (meas_package.timestamp_ - time_us_)/1000000.0;
      Prediction(dt);
    }
    time_us_ = meas_package.timestamp_;

	  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
      UpdateLidar(meas_package);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
      UpdateRadar(meas_package);
    }

}


void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  // Compute sigma points
  // create augmented mean vector
  std::cout << "prediction start\n";
  VectorXd x_aug = VectorXd(n_aug_);

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // create augmented mean state
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;
  // create augmented covariance matrix
  MatrixXd Q(2, 2);
  Q << std_a_*std_a_, 0,
       0, std_yawdd_*std_yawdd_;
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(n_aug_-n_x_, n_aug_-n_x_) = Q;
  //create square root matrix
  MatrixXd A(n_aug_,n_aug_);
  A = P_aug.llt().matrixL();
  // create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  float coeff = sqrt(lambda_aug_ + n_aug_);
  for(int i=1; i<n_aug_+1; i++){
      Xsig_aug.col(i) = x_aug + coeff*A.col(i-1);
  }
  for(int i=n_aug_+1; i<2*n_aug_+1; i++){
      Xsig_aug.col(i) = x_aug - coeff*A.col(i-n_aug_-1);
  }
  
  // predict sigma points
  for(int i=0; i<2*n_aug_+1; ++i){
    float px_k  = Xsig_aug(0,i);
    float py_k  = Xsig_aug(1,i);
    float v_k   = Xsig_aug(2,i);
    float phi_k = Xsig_aug(3,i);
    float phi_dot_k = Xsig_aug(4,i);
    float va_k = Xsig_aug(5,i);
    float vphidd_k = Xsig_aug(6,i);
    
    float p11;
    float p21;
    if(fabs(phi_dot_k) < 0.001){
        p11 = v_k*cos(phi_k)*delta_t;
        p21 = v_k*sin(phi_k)*delta_t;
    }else{
        p11 = (v_k/phi_dot_k)*(sin(phi_k + phi_dot_k*delta_t) - sin(phi_k));
        p21 = (v_k/phi_dot_k)*(-cos(phi_k + phi_dot_k*delta_t) + cos(phi_k));
    }
    float p12 = 0.5*delta_t*delta_t*cos(phi_k)*va_k;

    float p22 = 0.5*delta_t*delta_t*sin(phi_k)*va_k;
    float p31 = 0;
    float p32 = delta_t*va_k;
    float p41 = phi_dot_k*delta_t;
    float p42 = 0.5*delta_t*delta_t*vphidd_k;
    float p51 = 0;
    float p52 = delta_t*vphidd_k;
    
    float px_k_1    = px_k  + p11 + p12;
    float py_k_1    = py_k  + p21 + p22;
    float v_k_1     = v_k   + p31 + p32;
    float phi_k_1   = phi_k + p41 + p42; 
    float phi_dot_k_1 = phi_dot_k + p51 + p52;
    
    Xsig_pred_(0,i) = px_k_1;
    Xsig_pred_(1,i) = py_k_1;
    Xsig_pred_(2,i) = v_k_1;
    Xsig_pred_(3,i) = phi_k_1;
    Xsig_pred_(4,i) = phi_dot_k_1;
  }
  // predict state mea
  x_.fill(0.0);
  for(int i = 0; i < 2*n_aug_ + 1; ++i){
    x_ += weights_(i) * Xsig_pred_.col(i);  
  }

  // Predict covariance 
  P_.fill(0.0);
  for(int i = 0; i<2*n_aug_+1; ++i){
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    std::cout <<"i " << i << " x_diff: " << x_diff(3) << "\n";
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    P_ += weights_(i) * x_diff * x_diff.transpose() ;  
  }
  std::cout << "prediction start\n";
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  std::cout << "update lidar started\n";
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  int n_z_ = 2; // radar measurement dimension
  VectorXd z = VectorXd(n_z_);
  z <<
     meas_package.raw_measurements_[0],   // x
     meas_package.raw_measurements_[1];   // y
  MatrixXd H = MatrixXd(n_z_, n_x_); //
  H << 1,0,0,0,0,
       0,1,0,0,0;
       
  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  z_pred.fill(0.0);
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);
  MatrixXd R = MatrixXd(n_z_,n_z_);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;
  // predict measurement mean
  for(int i=0; i<2*n_aug_+1;++i){
    z_pred += H*Xsig_pred_.col(i) * weights_(i);  
  }
  // predict measurement cov
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = H*Xsig_pred_.col(i)-z_pred;

    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  S += R;
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = H * Xsig_pred_.col(i) - z_pred;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  // residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
  // angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  std::cout << "update lidar ended\n";
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  std::cout << "update radar start\n";
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  int n_z_ = 3; // radar measurement dimension
  VectorXd z = VectorXd(n_z_);
  z <<
     meas_package.raw_measurements_[0],   // rho in m
     meas_package.raw_measurements_[1],   // phi in rad
     meas_package.raw_measurements_[2];   // rho_dot in m/s

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);
  // compute
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;
    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                       // r
    Zsig(1,i) = atan2(p_y,p_x);                                // phi
    Zsig(2,i) = (p_x!=0 || p_y!=0)?(p_x*v1 + p_y*v2) / Zsig(0,i):0;   // r_dot
  }

  // mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; ++i) {
    z_pred = z_pred + Zsig.col(i) * weights_(i);
  }

  // innovation covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  MatrixXd R = MatrixXd(n_z_,n_z_);
  R <<  std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0,std_radrd_*std_radrd_;
  S = S + R;
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

  // angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  std::cout << "update radar ended\n";
}