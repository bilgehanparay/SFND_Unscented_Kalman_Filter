#ifndef UKF_H
#define UKF_H

#include "Eigen/Dense"
#include "measurement_package.h"

class UKF {
 public:
  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);


  // initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_; // todo

  // if this is false, laser measurements will be ignored (except for init)
  bool use_laser_; // done

  // if this is false, radar measurements will be ignored (except for init)
  bool use_radar_; // done

  // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  Eigen::VectorXd x_; // done

  // state covariance matrix
  Eigen::MatrixXd P_; // done

  // predicted sigma points matrix
  Eigen::MatrixXd Xsig_pred_; // to do

  // time when the state is true, in us
  long long time_us_; // todo

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_; // done

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_; // done

  // Laser measurement noise standard deviation position1 in m
  double std_laspx_; // done

  // Laser measurement noise standard deviation position2 in m
  double std_laspy_; // done

  // Radar measurement noise standard deviation radius in m
  double std_radr_; // done

  // Radar measurement noise standard deviation angle in rad
  double std_radphi_; // done

  // Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ; // done

  // Weights of sigma points
  Eigen::VectorXd weights_; // done

  // State dimension
  int n_x_; // done

  // Augmented state dimension
  int n_aug_; // done

  // Sigma point spreading parameter
  double lambda_; // done
  double lambda_aug_;
};

#endif  // UKF_H