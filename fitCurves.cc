#include <Eigen/Eigenvalues>
#include "calcGraphProps.h"
#include "fitCurve.h"

using namespace std;
using namespace Eigen;

double test(double x) {
  return x;
}

fxs fitCurves::fitFx(const MatrixXd &xData,  const MatrixXd &yData, vector<fx> toFit) {
  fxs f;
  f.push_back(test);
  return f;
} 

