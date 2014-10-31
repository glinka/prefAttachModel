#include <Eigen/Dense>
#include <Eigen/LU>
#include <iostream>
#include "calcGraphProps.h"
#include "fitCurves.h"

using namespace std;
using namespace Eigen;

vector<double> fitCurves::fitFx(const vector<double> &xData, const vector<double> &yData, fxs toFit) {
    //data points are (xData[i], yData[i])
    //"vectors" are here considered to be columnar
    int nFns = toFit.size();
    int nPts = xData.size();
    int i, j;
    vector<double> v;
    MatrixXd fxsEval(nPts, nFns);
    MatrixXd yDataCpy(nPts, 1);
    //populate fxsEval with the evaluations of the f(x)s at the given points
    for(i = 0; i < nPts; i++) {
	for(j = 0; j < nFns; j++) {
	    fxsEval(i, j) = (*toFit[j])(xData[i]);
	}
	yDataCpy(i) = yData[i];
    }
    MatrixXd coeffs(nFns, 1);
    coeffs = (((((fxsEval.transpose())*fxsEval).lu()).inverse())*(fxsEval.transpose()))*yDataCpy;
    return vector<double>(coeffs.data(), coeffs.data() + nFns);
    /**    
	used when returning fnlx, but need coeffs to project so return vector<double>
	function<double(double)> fittedFns = [=] (double x) {
	double eval = 0;
	int count;
	for(count = 0; count < nFns; count++) {
	    eval += toFit[count](x)*coeffs(count);
	}
	return eval;
    };
    return fittedFns;
    **/
}

vector<double> fitCurves::fitFx_constrained(const vector<double>& time_data, const vector<double> &index_data, const vector< vector<double> > &ydata, fxs to_fit_coeff_evo, fxs to_fit_degs, const double tproj, const double integral) {
  // n coefficient data points
  // k variables (number of polynomial coefficients used to fit each coefficient's time evolution)
  // p sub-variables (number of polynomial coefficients in function used to fit degree distribution)
  const int n = time_data.size();
  const int k = to_fit_coeff_evo.size();
  const int p = to_fit_degs.size();
  const int nindices = index_data.size();

  MatrixXd A_sub(n, k);
  VectorXd y(n*p);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < k; j++) {
      A_sub(i, j) = (to_fit_coeff_evo[j])(time_data[i]);
    }
  }

  for(int i = 0; i < p; i++) {
    for(int j = 0; j < n; j++) {
      y(i*n + j) = ydata[i][j];
    }
  }

  VectorXd t_sub(k);
  for(int i = 0; i < k; i++) {
    t_sub(i) = (to_fit_coeff_evo[i])(tproj);
  }

  VectorXd w = VectorXd::Zero(p);
  for(int i = 0; i < p; i++) {
    for(int j = 0; j < nindices; j++) {
      w(i) += (to_fit_degs[i])(index_data[j]);
    }
  }

  // now that we have the pieces, combine into one monstrous matrix-vector equation
  MatrixXd A(n*p, k*p);
  VectorXd t(k*p);
  for(int i = 0; i < p; i++) {
    A.block(i*n, i*k, n, k) = A_sub;
    t.segment(i*k, k) = w(i)*t_sub;
  }
  
  // int c = 2*(integral - t.dot((A.transpose()*A).lu().inverse()*A.transpose()*y))/t.dot((A.transpose()*A).lu().inverse()*t);
  // VectorXd x = (A.transpose()*A).lu().inverse()*(c*t/2 + A.transpose()*y);
  VectorXd coeffs = (A.transpose()*A).lu().inverse()*((2*(integral - t.dot((A.transpose()*A).lu().inverse()*A.transpose()*y))/t.dot((A.transpose()*A).lu().inverse()*t))*t/2 + A.transpose()*y);
  // return vector<double>(coeffs.data(), coeffs.data() + m);
  vector<double> new_coeffs(p, 0);
  for(int i = 0; i < p; i++) {
    for(int j = 0; j < k; j++) {
      new_coeffs[i] += (to_fit_coeff_evo[j])(tproj);
    }
  }
  return new_coeffs;
    
  // // numbers get too large? normalize a bit to help
  // w.array() /= constraint;
  // int normalized_constraint = 1;
  // // applies constraint w^T*coeffs = constraint
  // // int c = 2*(normalized_constraint - w.dot((fxs_evals.transpose()*fxs_evals).lu().inverse()*fxs_evals.transpose()*y))/w.dot((fxs_evals.transpose()*fxs_evals).lu().inverse()*w);
  // // VectorXd coeffs = (fxs_evals.transpose()*fxs_evals).lu().inverse()*(w*c/2 + fxs_evals.transpose()*y);
  // VectorXd coeffs = (fxs_evals.transpose()*fxs_evals).lu().inverse()*(w*(2*(normalized_constraint - w.dot((fxs_evals.transpose()*fxs_evals).lu().inverse()*fxs_evals.transpose()*y))/(w.dot((fxs_evals.transpose()*fxs_evals).lu().inverse()*w)))/2 + fxs_evals.transpose()*y);
  // return vector<double>(coeffs.data(), coeffs.data() + m);
}
