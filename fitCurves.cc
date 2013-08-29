#include <Eigen/Dense>
#include <Eigen/LU>
#include "calcGraphProps.h"
#include "fitCurves.h"

using namespace std;
using namespace Eigen;

fnlx fitCurves::fitFx(const vector<double> &xData, const vector<double> &yData, fxs toFit) {
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
	    fxsEval(i, j) = toFit[j](xData[i]);
	}
	yDataCpy(i) = yData[i];
    }
    MatrixXd coeffs(nFns, 1);
    coeffs = (((((fxsEval.transpose())*fxsEval).lu()).inverse())*(fxsEval.transpose()))*yDataCpy;
    function<double(double)> fittedFns = [=] (double x) {
	double eval = 0;
	int count;
	for(count = 0; count < nFns; count++) {
	    eval += toFit[count](x)*coeffs(count);
	}
	return eval;
    };
    return fittedFns;
}
