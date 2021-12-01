#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;


Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse = VectorXd(4);
  rmse << 0,0,0,0;

  unsigned int t = estimations.size();
  
  if(estimations.size() == 0 || ground_truth.size() == 0){
     cout << "Estimations or ground truth vector had size 0" << endl;
     return rmse;
  }

  if (t!= ground_truth.size()){

    return rmse;
  }
  
  for(unsigned int n =0; n < estimations.size(); n++){
      VectorXd residuals = estimations[n] - ground_truth[n];
      residuals = residuals.array() * residuals.array();
      rmse += residuals;
   }

  
  
   rmse =  rmse / t;
   rmse = rmse.array().sqrt();

   return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */

  MatrixXd Hj = MatrixXd::Zero(3,4);
  
   float px = x_state(0); // x pos
   float py = x_state(1); // y pos
   float vx = x_state(2); // x velocity
   float vy = x_state(3); // y velocity

   // calculate px^2 and py^2
   float px_2 = px * px;
   float py_2 = py * py;


   // calculate repetitive operations
   float sum_pxy = px_2 + py_2; // sum px_2, py_2
   float pow_pxy2 = pow(sum_pxy, 3/2); // power 3/2 of usm px_2, py_2
   float aux = vy*px - vx*py;

   if((fabs(aux)<  0.00001) || (px == 0 || py == 0)){
      cout << "Can't divide by 0" << endl;
      return Hj;
   }



   Hj << px/sqrt(sum_pxy), py/sqrt(sum_pxy), 0,0,
         -(py)/sum_pxy, px/sum_pxy, 0, 0,
         py*aux/pow_pxy2, px*aux/pow_pxy2, px/sqrt(sum_pxy), py/sqrt(sum_pxy);

   return Hj;
   
}


MatrixXd Tools::CalculateCovariant(const double dt, const double noise_ax, const double noise_ay){

   MatrixXd Q_(4,4);

   float dt_2 = dt * dt;
   float dt_3 = dt_2 * dt;
   float dt_4 = dt_3 * dt;

   Q_ = MatrixXd(4, 4);
   Q_ << (dt_4/4)*noise_ax, 0, ((dt_3)/2)*noise_ax, 0,
      0, ((dt_4)/4)*noise_ay, 0, ((dt_3)/2)*noise_ay,
      ((dt_3)/2)*noise_ax, 0, (dt_2)*noise_ax, 0,
      0, ((dt_3)/2)*noise_ay, 0, (dt_2)*noise_ay;

   return Q_;

}
