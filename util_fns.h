#ifndef UTIL_FNS_H_
#define UTIL_FNS_H_

#include <vector>
#include <string>
#include <iostream>

namespace utils {

  double average(const std::vector<double>& v);
  double get_median(const std::vector<double>& v);
  std::vector<double> get_squared_distances(const std::vector< std::vector<double> >& vectors);
  double l2_norm(const std::vector<double>& x1, const std::vector<double>& x2);

  template <typename T>
    void save_matrix(const std::vector< std::vector<T> >& A, std::ofstream& output_file, const std::string header="", const char delim=',');
  template <typename T>
    void save_vector(const std::vector<T>& v, std::ofstream& output_file, const std::string header="", const char delim='\n');

} 

#include "util_fns.tpp"

#endif
