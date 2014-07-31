#include "util_fns.h"

namespace utils {

  double average(const std::vector<double>& v) {
    double average = 0;
    int n = v.size();
    for(int i = 0; i < n; i++) {
      average += v[i];
    }
    return average/n;
  }

} 
