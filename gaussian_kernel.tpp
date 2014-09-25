#include <cmath>
#include "gaussian_kernel.h"

Gaussian_Kernel::Gaussian_Kernel(const double epsilon): epsilon_(epsilon) {}

double Gaussian_Kernel::kernel(const std::vector<double>& x1, const std::vector<double>& x2) {
  // if(x1.size() != x2.size()) {
  //   std::cout << "array dimensions do not match" << std::endl;
  //   exit(1);
  // }
  int n = x1.size();
  double norm = 0;
  for(int i = 0; i < n; i++) {
    norm += std::pow(x1[i] - x2[i], 2);
  }
  return std::exp(-norm/epsilon_);
}
