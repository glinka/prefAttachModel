#ifndef FIT_CURVES_H
#define FIT_CURVES_H
#include <vector>

typedef double (*fx)(double);
typedef double (*fxy)(double, double);
typedef std::vector<fx> fxs;
//typedef std::vector<std::vector<double>> dMatrix;

class fitCurves {
 public:
  fitCurves() {};
  ~fitCurves() {};
  fxs fitFx(const std::vector<double> &xData, const std::vector<double> &yData, fxs toFit);
 private:
};
  
#endif
