#ifndef FIT_CURVES_H
#define FIT_CURVES_H
#include <Eigen/Dense>
#include <vector>

typedef double (*fx)(double);
typedef double (*fxy)(double);
typedef std::vector<fx> fxs;

class fitCurves {
 public:
  fitCurves() {};
  ~fitCurves() {};
  fxs fitFx(const Eigen::MatrixXd &xData, const Eigen::MatrixXd &yData, fxs toFit);
 private:
};
  
#endif
