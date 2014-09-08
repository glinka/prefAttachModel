#ifndef NEWTON_H_
#define NEWTON_H_

#include <Eigen/Dense>
#include "newton_wrapper.h"

class GMRES;
class Linear_Solver;
class pamCPI;

class Newton {
 public:
  Newton(const double tol_abs, const double tol_rel, const int itermax);
  ~Newton() {}
  Eigen::VectorXd find_zero(Eigen::VectorXd (*F)(const Eigen::VectorXd&), Eigen::MatrixXd (*DF)(const Eigen::VectorXd&), const Eigen::VectorXd& x0, const Linear_Solver& ls) const;
  Eigen::VectorXd find_zero(Eigen::VectorXd (*F)(const Eigen::VectorXd&), const Eigen::VectorXd& x0, const double dx, const GMRES& ls) const;
  Eigen::VectorXd find_zero(Eigen::VectorXd (*F)(const Eigen::VectorXd&, pamCPI), const Eigen::VectorXd& x0, const double dx, const GMRES& ls, pamCPI model) const;
 private:
  const double tol_abs_;
  const double tol_rel_;
  const int itermax_;
};

#endif
    
    
    
