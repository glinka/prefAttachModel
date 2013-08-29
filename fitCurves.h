#ifndef FIT_CURVES_H
#define FIT_CURVES_H
#include <vector>
#include <functional>

typedef double (*fx)(double);
typedef double (*fxy)(double, double);
typedef std::vector<fx> fxs;
typedef std::function<double(double)> fnlx;
//typedef std::vector<std::vector<double>> dMatrix;

class fitCurves {
 public:
  fitCurves() {};
  ~fitCurves() {};
  static std::vector<double> fitFx(const std::vector<double> &xData, const std::vector<double> &yData, fxs toFit);
 private:
};
  
#endif
