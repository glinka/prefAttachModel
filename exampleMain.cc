#include <iostream>
#include "fitCurve.h"

using namespace std;

double test2(double x) {
  return x;
}

int main(int argc, char *argv[]) {
  fitCurves f;
  int n = 2;
  Eigen::MatrixXd A(n, n);
  fxs fs;
  fs.push_back(test2);
  fxs g = f.fitFx(A, A, fs);
  cout << (*(g[0]))(2.5) << endl;
}
