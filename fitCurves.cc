#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
#include "calcGraphProps.h"
#include "fitCurve.h"

using namespace std;
using namespace Eigen;

fxs fitCurves::fitFx(const dMatrix &xData, const dMatrix &yData, fxs toFit) {
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
      fxsEval(i, j) = (*fxs[j])(xData[i]);
    }
    yDataCpy(i) = yData[i];
  }
  MatrixXd coeffs(nFns, 1);
  LUDecomposition<MatrixXd> lu((fxsEval.transpose())*fxsEval);
  coeffs = ((lu.inverse())*(fxsEval.transpose()))*yDataCpy;
  vector<double> coeffsVect;
  for(i = 0; i < nFns; i++) {
    coeffsVect.push_back(coeffs(i));
  }
  return f;
} 

