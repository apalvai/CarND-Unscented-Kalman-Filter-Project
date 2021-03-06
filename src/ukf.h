#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:
    
    ///* initially set to false, set to true in first call of ProcessMeasurement
    bool is_initialized_;
    
    ///* if this is false, laser measurements will be ignored (except for init)
    bool use_laser_;
    
    ///* if this is false, radar measurements will be ignored (except for init)
    bool use_radar_;
    
    ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
    VectorXd x_;
    
    ///* state covariance matrix
    MatrixXd P_;
    
    ///* predicted sigma points matrix
    MatrixXd Xsig_pred_;
    
    ///* time when the state is true, in us
    long long time_us_;
    
    ///* Process noise standard deviation longitudinal acceleration in m/s^2
    double std_a_;
    
    ///* Process noise standard deviation yaw acceleration in rad/s^2
    double std_yawdd_;
    
    ///* Laser measurement noise standard deviation position1 in m
    double std_laspx_;
    
    ///* Laser measurement noise standard deviation position2 in m
    double std_laspy_;
    
    ///* Radar measurement noise standard deviation radius in m
    double std_radr_;
    
    ///* Radar measurement noise standard deviation angle in rad
    double std_radphi_;
    
    ///* Radar measurement noise standard deviation radius change in m/s
    double std_radrd_ ;
    
    ///* Weights of sigma points
    VectorXd weights_;
    
    ///* State dimension
    int n_x_;
    
    ///* Augmented state dimension
    int n_aug_;
    
    ///* Sigma point spreading parameter
    double lambda_;
    
    
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
    
private:
    // previous timestamp
    long long previous_timestamp_;
    
    // process covariance matrix
    MatrixXd Q_;
    
    // augmented matrix
    MatrixXd Xsig_aug_;
    
    // Sigma point spreading parameter for augmented
    double lambda_aug_;
    
    int n_sigma_;
    
    // lidar measurement covariance matrix
    MatrixXd R_lidr_;
    
    // radar measurement covariance matrix
    MatrixXd R_radr_;
    
    // tools
    Tools tools;
    
    // NIS_lidr_
    double NIS_lidr_;
    
    // NIS_radr_
    double NIS_radr_;
    
    /**
     * Used to generate Sigma points including process noise
     */
    void AugmentedSigmaPoints(MatrixXd* Xsig_out);
    
    /**
     * Used to apply CTRV transformation for predicting sigma points through state transition
     */
    void ApplyCTRVTransform(double delta_t);
    
    /**
     * Used to compute state mean and coveraince based on predicted sigma points
     */
    void PredictStateMeanAndCovariance();
    
    /**
     * Used to transform predicted sigma points to RADAR measurement space
     */
    void CreateRadarMeasurementSigmaPoints(MatrixXd* Zsig_out);
    
    /**
     * Used to transform predicted sigma points to LASER measurement space
     */
    void CreateLaserMeasurementSigmaPoints(MatrixXd* Zsig_out);
};

#endif /* UKF_H */
