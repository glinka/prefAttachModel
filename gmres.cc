#include <cmath>
// #include "la_ops.h"
#include <Eigen/Dense>
#include "gmres.h"
#include "iters_exception.h"
#include "pamCPI.h"

// TESTING
#include <iostream>

GMRES::GMRES(const double tol, const double kmax): tol_(tol), kmax_(kmax) {}

Eigen::VectorXd GMRES::solve_linear_system(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const Eigen::VectorXd& x0) const {
  const int n = b.size();
  Eigen::VectorXd r = b - A*x0;
  /*
    V will store the basis vectors for
    the Krylov subspace as column vectors
  */
  Eigen::MatrixXd V(n, kmax_ + 1);
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(kmax_ + 1, kmax_);
  // V[0] is really x0, and not part of the Krylov subspace
  V.col(0) = r/r.norm();
  double rho = r.norm();
  const double b_norm = b.norm();
  int k = -1;
  Eigen::VectorXd g(kmax_ + 1);
  Eigen::VectorXd h(kmax_ + 1);
  g(0) = rho;
  double c[kmax_];
  double s[kmax_];
  Eigen::MatrixXd G = Eigen::MatrixXd::Identity(kmax_ + 1, kmax_ + 1);
  while(rho > tol_*b_norm && k < kmax_-1) {
    k++;
    V.col(k+1) = A*V.col(k);
    for(int j = 0; j < k+1; j++) {
      h(j) = V.col(k+1).dot(V.col(j));
      V.col(k+1) -= h(j)*V.col(j);
    }
    h(k+1) = V.col(k+1).norm();
    // test for orthogonality
    V.col(k+1) = V.col(k+1)/h(k+1);
    // perform first loop separately to avoid if?
    // will not execute during first loop, as desired
    for(int i = 0; i < k; i++) {
      double temp1 = c[i]*h(i) - s[i+1]*h(i+1);
      double temp2 = s[i]*h(i) + c[i+1]*h(i+1);
      h(i) = temp1;
      h(i+1) = temp2;
      // executed every time, wasteful
      G(k-1, k-1) = 1;
      G(k-1, k) = 0;
      G(k, k-1) = 0;
    }
    double nu = sqrt(h(k)*h(k) + h(k+1)*h(k+1));
    c[k] = h(k)/nu;
    s[k] = -h(k+1)/nu;
    h(k) = c[k]*h(k) - s[k]*h(k+1);
    h(k+1) = 0;
    // set Givens matrix
    G(k, k) = c[k];
    G(k, k+1) = -s[k];
    G(k+1, k) = s[k];
    G(k+1, k+1) = c[k];
    g = G*g;
    for(int i = 0; i < k+1; i++) {
      H(i, k) = h(i);
    }
    rho = abs(g(k+1));
  }

  Iters_Exception::test_iters(k, kmax_ - 1);

  Eigen::VectorXd y(k+1);
  // y(k+1) = 1;
  // solve for y^k
  for(int i = k; i >= 0; i--) {
    double sum = 0;
    for(int j = k; j > i; j--) {
      sum += y(j)*H(i,j);
    }
    y(i) = (g(i) - sum)/H(i,i);
  }

  Eigen::VectorXd x = x0 + V.leftCols(k+1)*y;

  return x;

  // for(int i = 0; i < n; i++) {
  //   double sum = 0;
  //   for(int j = 0; j < k+2; j++) {
  //     // because V is transposed
  //     sum += V(j,i)*y(j);
  //   }
  //   x(i) = x0(i) + sum;
  // }

  // cleanup, can move to one line?
  // delete r;
  // delete V;
  // delete H;
  // delete g;
  // delete h;
  // delete y;
}

namespace la {

  Eigen::VectorXd v_next(const Eigen::VectorXd& x, const Eigen::VectorXd& w, const double h, Eigen::VectorXd (*F)(const Eigen::VectorXd&)) {
    return w.norm()*(F(x + h*x.norm()*w/(w.norm())) - F(x))/(h*x.norm());
  }

  // for use when x = 0
  Eigen::VectorXd v_next_init(const Eigen::VectorXd& x, const Eigen::VectorXd& w, const double h, Eigen::VectorXd (*F)(const Eigen::VectorXd&)) {
    return w.norm()*(F(h*w/(w.norm())) - F(x))/h;
  }
}


Eigen::VectorXd GMRES::solve_linear_system(Eigen::VectorXd (*F)(const Eigen::VectorXd&), const Eigen::VectorXd& x, const Eigen::VectorXd& x0, const double dx) const {
  Eigen::VectorXd b = F(x);
  const int n = b.size();
  /*
    V will store the basis vectors for
    the Krylov subspace as column vectors
  */
  Eigen::MatrixXd V(n, kmax_ + 1);
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(kmax_ + 1, kmax_);
  // V[0] is really x0, and not part of the Krylov subspace
  Eigen::VectorXd (*v_next)(const Eigen::VectorXd& x, const Eigen::VectorXd& w, const double h, Eigen::VectorXd (*F)(const Eigen::VectorXd&));
  Eigen::VectorXd r;
  if(x0.norm() == 0) {
    v_next = la::v_next_init;
    r = b;
  }
  else {
    v_next = la::v_next;
    r = b - v_next(x0, x0, dx, F);
  }
  V.col(0) = r/r.norm();
  double rho = r.norm();
  const double b_norm = b.norm();
  int k = -1;
  Eigen::VectorXd g(kmax_ + 1);
  Eigen::VectorXd h(kmax_ + 1);
  g(0) = rho;
  double c[kmax_];
  double s[kmax_];
  Eigen::MatrixXd G = Eigen::MatrixXd::Identity(kmax_ + 1, kmax_ + 1);
  // select appropriate method for determining v_{k+1}
  while(rho > tol_*b_norm && k < kmax_-1) {
    k++;
    V.col(k+1) = v_next(x, V.col(k), dx, F);
    for(int j = 0; j < k+1; j++) {
      h(j) = V.col(k+1).dot(V.col(j));
      V.col(k+1) -= h(j)*V.col(j);
    }
    h(k+1) = V.col(k+1).norm();
    // test for orthogonality
    V.col(k+1) = V.col(k+1)/h(k);
    // perform first loop separately to avoid if?
    // will not execute during first loop, as desired
    for(int i = 0; i < k; i++) {
      double temp1 = c[i]*h(i) - s[i+1]*h(i+1);
      double temp2 = s[i]*h(i) + c[i+1]*h(i+1);
      h(i) = temp1;
      h(i+1) = temp2;
      // executed every time, wasteful
      G(k-1, k-1) = 1;
      G(k-1, k) = 0;
      G(k, k-1) = 0;
    }
    double nu = sqrt(h(k)*h(k) + h(k+1)*h(k+1));
    c[k] = h(k)/nu;
    s[k] = -h(k+1)/nu;
    h(k) = c[k]*h(k) - s[k]*h(k+1);
    h(k+1) = 0;
    // set Givens matrix
    G(k, k) = c[k];
    G(k, k+1) = -s[k];
    G(k+1, k) = s[k];
    G(k+1, k+1) = c[k];
    g = G*g;
    for(int i = 0; i < k+1; i++) {
      H(i, k) = h(i);
    }
    rho = abs(g(k+1));
  }
  Eigen::VectorXd y(k+2);
  y(k+1) = 1;
  // solve for y^k
  for(int i = k; i >= 0; i--) {
    double sum = 0;
    for(int j = k; j > i; j--) {
      sum += y(j)*H(i,j);
    }
    y(i) = (g(i) - sum)/H(i,i);
  }

  return x0 + V.leftCols(k+2)*y;

  // for(int i = 0; i < n; i++) {
  //   double sum = 0;
  //   for(int j = 0; j < k+2; j++) {
  //     // because V is transposed
  //     sum += V(j,i)*y(j);
  //   }
  //   x(i) = x0(i) + sum;
  // }

  // cleanup, can move to one line?
  // delete r;
  // delete V;
  // delete H;
  // delete g;
  // delete h;
  // delete y;
}

// FROM HERE ON
// are functions specifically for pamCPI

namespace la {

  Eigen::VectorXd v_next(const Eigen::VectorXd& x, const Eigen::VectorXd& w, const double h, Eigen::VectorXd (*F)(const Eigen::VectorXd&, pamCPI), pamCPI model) {
    return w.norm()*(F(x + h*x.norm()*w/(w.norm()), model) - F(x, model))/(h*x.norm());
  }

  // for use when x = 0

  Eigen::VectorXd v_next_zero(const Eigen::VectorXd& x, const Eigen::VectorXd& w, const double h, Eigen::VectorXd (*F)(const Eigen::VectorXd&, pamCPI), pamCPI model) {
    return w.norm()*(F(h*w/(w.norm()), model) - F(x, model))/h;
  }

  // "init" because often we will be solving DF(x) dx = F(x) where
  // the initial guess for dx is 0
  Eigen::VectorXd v_next_init(const Eigen::VectorXd& x, const Eigen::VectorXd& w, const double h, Eigen::VectorXd (*F)(const Eigen::VectorXd&, pamCPI), pamCPI model) {
    return Eigen::VectorXd::Zero(x.size());
  }
}

Eigen::VectorXd GMRES::solve_linear_system(Eigen::VectorXd (*F)(const Eigen::VectorXd&, pamCPI), const Eigen::VectorXd& x, const Eigen::VectorXd& x0, const double dx, pamCPI model) const {
  Eigen::VectorXd b = F(x, model);
  const int n = b.size();
  /*
    V will store the basis vectors for
    the Krylov subspace as column vectors
  */
  Eigen::MatrixXd V(n, kmax_ + 1);
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(kmax_ + 1, kmax_);
  // V[0] is really x0, and not part of the Krylov subspace
  Eigen::VectorXd (*v_next)(const Eigen::VectorXd& x, const Eigen::VectorXd& w, const double h, Eigen::VectorXd (*F)(const Eigen::VectorXd&, pamCPI), pamCPI model);
  Eigen::VectorXd r;
  // assume w.norm (= V.col(k).norm) \neq 0
  // otherwise we should have ended iterations
  // earlier
  if (x.norm() == 0) {
    v_next = la::v_next_zero;
  }
  else {
    v_next = la::v_next;
  }
  // however, x0.norm could still be zero, which
  // we account for below
  if (x0.norm() == 0) {
    r = b;
  }
  else {
    r = b - v_next(x, x0, dx, F, model);
  }    
  V.col(0) = r/r.norm();
  std::cout << "V_0 sum: " << V.col(0).sum() << std::endl;
  double rho = r.norm();
  const double b_norm = b.norm();
  int k = -1;
  Eigen::VectorXd g(kmax_ + 1);
  Eigen::VectorXd h(kmax_ + 1);
  g(0) = rho;
  double c[kmax_];
  double s[kmax_];
  Eigen::MatrixXd G = Eigen::MatrixXd::Identity(kmax_ + 1, kmax_ + 1);
  // select appropriate method for determining v_{k+1}
  while(rho > tol_*b_norm && k < kmax_-1) {
    k++;
    V.col(k+1) = v_next(x, V.col(k), dx, F, model);
    std::cout << "V_" << k+2 << " 1-norm: " << V.col(k+1).sum() << std::endl;
    for(int j = 0; j < k+1; j++) {
      h(j) = V.col(k+1).dot(V.col(j));
      V.col(k+1) -= h(j)*V.col(j);
    }
    h(k+1) = V.col(k+1).norm();
    // test for orthogonality
    V.col(k+1) = V.col(k+1)/h(k);
    // perform first loop separately to avoid if?
    // will not execute during first loop, as desired
    for(int i = 0; i < k; i++) {
      double temp1 = c[i]*h(i) - s[i+1]*h(i+1);
      double temp2 = s[i]*h(i) + c[i+1]*h(i+1);
      h(i) = temp1;
      h(i+1) = temp2;
      // executed every time, wasteful
      G(k-1, k-1) = 1;
      G(k-1, k) = 0;
      G(k, k-1) = 0;
    }
    double nu = sqrt(h(k)*h(k) + h(k+1)*h(k+1));
    c[k] = h(k)/nu;
    s[k] = -h(k+1)/nu;
    h(k) = c[k]*h(k) - s[k]*h(k+1);
    h(k+1) = 0;
    // set Givens matrix
    G(k, k) = c[k];
    G(k, k+1) = -s[k];
    G(k+1, k) = s[k];
    G(k+1, k+1) = c[k];
    g = G*g;
    for(int i = 0; i < k+1; i++) {
      H(i, k) = h(i);
    }
    rho = abs(g(k+1));
  }
  Eigen::VectorXd y(k+2);
  y(k+1) = 1;
  // solve for y^k
  for(int i = k; i >= 0; i--) {
    double sum = 0;
    for(int j = k; j > i; j--) {
      sum += y(j)*H(i,j);
    }
    y(i) = (g(i) - sum)/H(i,i);
  }

  return x0 + V.leftCols(k+2)*y;
}
