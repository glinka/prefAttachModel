#include <Eigen/Dense>
#include "newton.h"
#include "gmres.h"
#include "iters_exception.h"
#include "pamCPI.h"

// TESTING
#include <iostream>


Newton::Newton(const double tol_abs, const double tol_rel, const int itermax): tol_abs_(tol_abs), tol_rel_(tol_rel), itermax_(itermax) {}

Eigen::VectorXd Newton::find_zero(Eigen::VectorXd (*F)(const Eigen::VectorXd&), Eigen::MatrixXd (*DF)(const Eigen::VectorXd&), const Eigen::VectorXd& x0, const Linear_Solver& ls) const {
  const int n = x0.size();
  double r0 = F(x0).norm();
  double r = r0;
  Eigen::VectorXd zeros = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd x = x0;
  int iters = 0;
  while(r > r0*tol_rel_ + tol_abs_ && iters < itermax_) {
    // default to initial x = {0, 0, ..., 0}
    x -= ls.solve_linear_system(DF(x), F(x), zeros);
    r = F(x).norm();
    iters++;
  }
  // probably a better way to do this
  Iters_Exception::test_iters(iters, itermax_);
  return x;
}

Eigen::VectorXd Newton::find_zero(Eigen::VectorXd (*F)(const Eigen::VectorXd&), const Eigen::VectorXd& x0, const double dx, const GMRES& ls) const {
  const int n = x0.size();
  double r0 = F(x0).norm();
  double r = r0;
  Eigen::VectorXd zeros = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd x = x0;
  int iters = 0;
  while(r > r0*tol_rel_ + tol_abs_ && iters < itermax_) {
    // default to initial x = {0, 0, ..., 0}
    x += ls.solve_linear_system(F, x, zeros, dx);
    r = F(x).norm();
    iters++;
  }
  // probably a better way to do this
  Iters_Exception::test_iters(iters, itermax_);
  return x;
}

Eigen::VectorXd Newton::find_zero(Eigen::VectorXd (*F)(const Eigen::VectorXd&, pamCPI), const Eigen::VectorXd& x0, const double dx, const GMRES& ls, pamCPI model) const {
  const int n = x0.size();
  double r0 = F(x0, model).norm();
  double r = r0;
  Eigen::VectorXd zeros = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd x = x0;
  int iters = 0;
  while(r > r0*tol_rel_ + tol_abs_ && iters < itermax_) {
    // default to initial x = {0, 0, ..., 0}
    std::cout << "newton residual: " << r << std::endl;
    Eigen::VectorXd dx_newton = ls.solve_linear_system(F, x, zeros, dx, model);
    // do a line search
    while((x+dx_newton).minCoeff() < 0) {
      dx_newton /= 2;
    }
    x += dx_newton;
    r = F(x, model).norm();
    iters++;
  }
  // probably a better way to do this
  Iters_Exception::test_iters(iters, itermax_);
  return x;
}
