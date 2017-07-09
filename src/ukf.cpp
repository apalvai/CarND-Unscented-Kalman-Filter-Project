#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

#include "measurement_package.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

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
    P_ = MatrixXd(5, 5);
    
    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 3;
    
    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 3;
    
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
     TODO:
     
     Complete the initialization. See ukf.h for other member properties.
     
     Hint: one or more values initialized above might be wildly off...
     */
    
    n_x_ = 5;
    n_aug_ = 7;
    n_sigma_ = 2*n_aug_+1;
    lambda_ = 3 - n_x_;
    lambda_aug_ = 3 - n_aug_;
    
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0 ,0, 0, 1;
    
    Q_ = MatrixXd(2, 2);
    Q_ << std_a_*std_a_, 0,
          0, fmod(std_yawdd_*std_yawdd_, 2*M_PI);
    cout << "Q_: " << endl << Q_ <<endl;
    
    Xsig_aug_ = MatrixXd(n_aug_, n_sigma_);
    Xsig_aug_.fill(0);
    
    //create matrix with predicted sigma points as columns
    Xsig_pred_ = MatrixXd(n_x_, n_sigma_);
    Xsig_pred_.fill(0);
    
    //create vector for weights
    weights_ = VectorXd(n_sigma_);
    
    double weight_1 = lambda_aug_/(lambda_aug_+n_aug_);
    double weight_n = 1/((lambda_aug_+n_aug_)*2);
    
    //set weights
    for (int i=0; i<n_sigma_; i++)
    {
        if(i == 0) {
            weights_(i) = weight_1;
        } else {
            weights_(i) = weight_n;
        }
    }
    cout << "weights: " << endl << weights_ <<endl;
    
    // lidar measurement covariance matrix
    R_lidr_ = MatrixXd(2, 2);
    R_lidr_ << std_laspx_*std_laspx_, 0,
               0, std_laspy_*std_laspy_;
    cout << "R_lidr_: " << endl << R_lidr_ <<endl;
    
    // radar measurement covariance matrix
    R_radr_ = MatrixXd(3, 3);
    R_radr_ << std_radr_*std_radr_, 0, 0,
               0, fmod(std_radphi_*std_radphi_, 2*M_PI), 0,
               0, 0, std_radrd_*std_radrd_;
    cout << "R_radr_: " << endl << R_radr_ <<endl;
    
    // Create an instance of tools
    tools = Tools();
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    
    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    
    if (is_initialized_ == false) {
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            /**
             Convert radar measurement into UKF state vector and initialize it.
             */
            VectorXd radarMeasurementVector = VectorXd(3);
            radarMeasurementVector << meas_package.raw_measurements_[0],
                                      meas_package.raw_measurements_[1],
                                      meas_package.raw_measurements_[2];
            x_ << tools.ConvertRadarMeasurementToUKFStateVector(radarMeasurementVector);
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            /**
             Initialize state.
             */
            x_ << meas_package.raw_measurements_[0],
                  meas_package.raw_measurements_[1],
                  0,
                  0,
                  0;
        }
        
        previous_timestamp_ = meas_package.timestamp_;
        
        // done initializing, no need to predict or update
        is_initialized_ = true;
        
        cout << "initialized UKF with state: " << endl << x_ << endl;
        return;
    }
    
    if (use_laser_ == false && meas_package.sensor_type_ == MeasurementPackage::LASER) {
        cout << "Skip processing laser measurements." << endl;
        return;
    } else if (use_radar_ == false && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        cout << "Skip processing radar measurements." << endl;
        return;
    }
    
    /*****************************************************************************
     *  Sigma Point Generation
     ****************************************************************************/
    
    AugmentedSigmaPoints(&Xsig_aug_);
    cout << "created sigma points: " << endl << Xsig_aug_ << endl;
    
    /*****************************************************************************
     *  Prediction
     ****************************************************************************/
    
    double delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0; //time diff in sec
    previous_timestamp_ = meas_package.timestamp_;
    Prediction(delta_t);
    cout << "predicted state: " << endl << x_ << endl;
    cout << "predicted covariance: " << endl << P_ << endl;
    cout << "predicted sigma points: " << endl << Xsig_pred_ << endl;
    
    /*****************************************************************************
     *  Update
     ****************************************************************************/
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        UpdateRadar(meas_package);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
        UpdateLidar(meas_package);
    }
    
    cout << "updated state after measurement: " << endl << x_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    /**
     Estimate the object's location. Modify the state
     vector, x_. Predict sigma points, the state, and the state covariance matrix.
     */
    
    ApplyCTRVTransform(delta_t);
    
    PredictStateMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    /**
     TODO:
     
     Complete this function! Use lidar data to update the belief about the object's
     position. Modify the state vector, x_, and covariance, P_.
     
     You'll also need to calculate the lidar NIS.
     */
    
    int n_z = 2;
    
    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, n_sigma_);
    CreateLaserMeasurementSigmaPoints(&Zsig);
    
    //calculate mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred = Zsig*weights_;
    cout << "laser z_pred: " << endl << z_pred << endl;
    
    //calculate measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0);
    for (int i=0; i<n_sigma_; i++)
    {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }
    S = S + R_lidr_;
    
    // measurement received
    VectorXd z = VectorXd(n_z);
    z << meas_package.raw_measurements_[0],
         meas_package.raw_measurements_[1];
    cout << "laser measurement point: " << endl << z << endl;
    
    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0);
    for(int i=0; i<n_sigma_; i++)
    {
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        
        VectorXd z_diff = Zsig.col(i) - z_pred;
        
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
    
    //calculate Kalman gain K
    MatrixXd K = MatrixXd(n_x_, n_z);
    K = Tc * S.inverse();
    // std::cout << "K: " << std::endl << K << std::endl;
    
    //update state mean and covariance matrix
    x_ = x_ + K * (z - z_pred);
    
    P_ = P_ - K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    /**
     Use radar data to update the belief about the object's
     position. Modify the state vector, x_, and covariance, P_.
     
     You'll also need to calculate the radar NIS.
     */
    
    int n_z = 3;
    
    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, n_sigma_);
    CreateRadarMeasurementSigmaPoints(&Zsig);
    
    //calculate mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred = Zsig*weights_;
    cout << "radar z_pred: " << endl << z_pred << endl;
    
    //calculate measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0);
    for (int i=0; i<n_sigma_; i++)
    {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        z_diff(1) = fmod(z_diff(1), 2*M_PI);
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }
    
    S = S + R_radr_;
    
    // measurement received
    VectorXd z = VectorXd(n_z);
    z << meas_package.raw_measurements_[0],
         meas_package.raw_measurements_[1],
         meas_package.raw_measurements_[2];
    cout << "radar measurement point: " << endl << z << endl;
    
    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0);
    for(int i=0; i<n_sigma_; i++)
    {
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        x_diff(3) = fmod(x_diff(3), 2*M_PI);
        
        VectorXd z_diff = Zsig.col(i) - z_pred;
        z_diff(1) = fmod(z_diff(1), 2*M_PI);
        
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
    
    //calculate Kalman gain K
    MatrixXd K = MatrixXd(n_x_, n_z);
    K = Tc * S.inverse();
    // std::cout << "K: " << std::endl << K << std::endl;
    
    //update state mean and covariance matrix
    x_ = x_ + K * (z - z_pred);
    x_(3) = fmod(x_(3), 2*M_PI);
    
    P_ = P_ - K * S * K.transpose();
}


#pragma mark - Private functions

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {
    
    //create augmented mean state
    VectorXd x_aug = VectorXd(n_aug_);
    x_aug << x_,
             0,
             0;
    cout << "x_aug: " << endl << x_aug << endl;
    
    //create augmented state covariance
    MatrixXd zero5_2 = MatrixXd(5, 2);
    zero5_2.fill(0);
    
    MatrixXd zero2_5 = MatrixXd(2, 5);
    zero2_5.fill(0);
    
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    P_aug << P_, zero5_2,
             zero2_5, Q_;
    
    cout << "P_aug: " << endl << P_aug << endl;
    
    //create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma_);
    Xsig_aug.fill(0);
    
    //calculate square root of P
    MatrixXd A_aug = P_aug.llt().matrixL();
    
    // compute sqrt((lambda_aug_+n_aug)*P)
    MatrixXd diff = sqrt(lambda_aug_ + n_aug_) * A_aug;
    cout << "diff: " << endl << diff << endl;
    
    for (int i=0; i<n_sigma_; i++)
    {
        VectorXd col = VectorXd(n_aug_);
        col.fill(0);
        if (i == 0) {
            col = x_aug;
        } else if (i < n_aug_ + 1) {
            // cout << "diff.col(i-1): " << endl << diff.col(i-1) << endl;
            col = x_aug + diff.col(i-1);
        } else {
            // cout << "diff.col(i-n_aug_-1): " << endl << diff.col(i-n_aug_-1) << endl;
            col = x_aug - diff.col(i-n_aug_-1);
        }
        
        col(3) = fmod(col(3), 2*M_PI);
        
        Xsig_aug.col(i) = col;
    }
    
    *Xsig_out = Xsig_aug;
}

void UKF::ApplyCTRVTransform(double delta_t) {
    
    double delta_t2_by_2 = delta_t*delta_t/2;
    
    //predict sigma points
    for (int i=0; i<2*n_aug_+1; i++)
    {
        VectorXd col = VectorXd(n_aug_);
        col << Xsig_aug_.col(i);
        
        VectorXd x_k = VectorXd(n_x_);
        x_k << col(0),
               col(1),
               col(2),
               col(3),
               col(4);
        
        double cos_si = cos(x_k(3));
        double sin_si = sin(x_k(3));
        
        VectorXd x_k_step = VectorXd(n_x_);
        if (x_k(4) == 0)
        {
            x_k_step << x_k(2)*cos_si*delta_t,
                        x_k(2)*sin_si*delta_t,
                        0,
                        0,
                        0;
        }
        else
        {
            double v_si_k = x_k(2)/x_k(4);
            double si_dot_dt = x_k(4)*delta_t;
            double cos_ssdot = cos(x_k(3)+si_dot_dt);
            double sin_ssdot = sin(x_k(3)+si_dot_dt);
            
            x_k_step << v_si_k*(sin_ssdot - sin_si),
                        v_si_k*(cos_si - cos_ssdot),
                        0,
                        si_dot_dt,
                        0;
        }
        
        VectorXd x_k_noise = VectorXd(n_x_);
        x_k_noise << delta_t2_by_2*cos_si*col(5),
                     delta_t2_by_2*sin_si*col(5),
                     delta_t*col(5),
                     delta_t2_by_2*col(6),
                     delta_t*col(6);
        
        Xsig_pred_.col(i) = x_k + x_k_step + x_k_noise;
        
        Xsig_pred_.col(i)(3) = fmod(Xsig_pred_.col(i)(3), 2*M_PI);
    }
}

void UKF::PredictStateMeanAndCovariance() {
    
    //compute mean of the predicted sigma points
    x_ = Xsig_pred_ * weights_;
    x_(3) = fmod(x_(3), 2*M_PI);
    
    //compute covariance of the predicted sigma points
    MatrixXd X_diff = MatrixXd(n_x_, n_sigma_);
    MatrixXd X_diff_w = MatrixXd(n_x_, n_sigma_);
    
    for (int i=0; i<n_sigma_; i++)
    {
        X_diff.col(i) = Xsig_pred_.col(i) - x_;
        
        //angle normalization
        X_diff.col(i)(3) = fmod(X_diff.col(i)(3), 2*M_PI);
        
        X_diff_w.col(i) = weights_(i)*X_diff.col(i);
        X_diff_w.col(i)(3) = fmod(X_diff_w.col(i)(3), 2*M_PI);
    }
    
    P_ = X_diff_w * X_diff.transpose();
}

void UKF::CreateRadarMeasurementSigmaPoints(MatrixXd* Zsig_out) {
    
    int n_z = 3;
    
    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, n_sigma_);
    Zsig.fill(0);
    
    //transform sigma points into measurement space
    for (int i=0; i<n_sigma_; i++)
    {
        VectorXd col = Xsig_pred_.col(i);
        
        double px = col(0);
        double py = col(1);
        double v  = col(2);
        double si = col(3);
        double c1 = sqrt(px*px + py*py);
        
        if (abs(px) < 0.0001)
        {
            std::cout << "ERROR:" << "px: " << px << std::endl;
            continue;
        }
        
        if (abs(c1) < 0.0001)
        {
            std::cout << "ERROR: " << "px: " << px << "py: " << py << std::endl;
            continue;
        }
        
        double angle = atan2(py, px);
        angle = fmod(angle, 2*M_PI);
        
        double c2 = (px*cos(si)+py*sin(si))*v/c1;
        
        Zsig.col(i) << c1,
        angle,
        c2;
    }
    cout << "predicted measurement sigma points (radar): " << endl << Zsig << endl;
    *Zsig_out = Zsig;
}

void UKF::CreateLaserMeasurementSigmaPoints(MatrixXd* Zsig_out) {
    
    int n_z = 2;
    
    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, n_sigma_);
    Zsig.fill(0);
    
    //transform sigma points into measurement space
    for (int i=0; i<n_sigma_; i++)
    {
        double px = Xsig_pred_.col(i)(0);
        double py = Xsig_pred_.col(i)(1);
        
        Zsig.col(i) << px,
                       py;
    }
    cout << "predicted measurement sigma points (laser): " << endl << Zsig << endl;
    *Zsig_out = Zsig;
}
