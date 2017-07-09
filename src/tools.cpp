#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    rmse << 0,0,0,0;
    
    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    if(estimations.size() != ground_truth.size()
       || estimations.size() == 0){
        cout << "Invalid estimation or ground_truth data" << endl;
        return rmse;
    }
    
    //accumulate squared residuals
    for(unsigned int i=0; i < estimations.size(); ++i) {
        
        VectorXd residual = estimations[i] - ground_truth[i];
        
        //coefficient-wise multiplication
        residual = residual.array()*residual.array();
        rmse += residual;
    }
    
    //calculate the mean
    rmse = rmse/estimations.size();
    
    //calculate the squared root
    rmse = rmse.array().sqrt();
    
    cout << "rmse: " << endl << rmse << endl;
    
    //return the result
    return rmse;
}

VectorXd Tools::ConvertRadarMeasurementToUKFStateVector(const VectorXd& x_radar) {
    
    VectorXd x_ukf(5);
    x_ukf.fill(0);
    
    std::cout << "x_radar: " << std::endl << x_radar << std::endl;
    
    //recover state parameters
    float range = x_radar(0);
    float bearing = x_radar(1);
    float range_rate = x_radar(2);
    
    //pre-compute a set of terms to avoid repeated calculation
    float px = range * cos(bearing);
    float py = range * sin(bearing);
    
    x_ukf << px,
             py,
             range_rate,
             bearing,
             0;
    
    return x_ukf;
}
