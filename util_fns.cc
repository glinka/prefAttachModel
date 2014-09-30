#include <algorithm>
#include <cmath>
#include <sstream>
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

  double get_median(const std::vector<double>& v) {
    std::vector<double> vcopy = v;
    std::sort(vcopy.begin(), vcopy.end());
    int n = vcopy.size();
    if(n % 2 == 0) {
      return (vcopy[n/2] + vcopy[n/2 + 1])/2.0;
    }
    else {
      return vcopy[n/2 + 1];
    }
  }

  std::vector<double> get_squared_distances(const std::vector< std::vector<double> >& vectors) {
    const int nvects = vectors.size();
    int ncombos = (nvects*(nvects-1))/2;
    std::vector<double> squared_distances(ncombos);
    int counter = 0;
    for(int i = 0; i < nvects; i++) {
      for(int j = i+1; j < nvects; j++) {
	squared_distances[counter++] = pow(l2_norm(vectors[i], vectors[j]), 2);
      }
    }
    return squared_distances;
  }
  
  double l2_norm(const std::vector<double>& x1, const std::vector<double>& x2) {
    double norm = 0;
    for(int i = 0; i < x1.size(); i++) {
      norm += pow(x1[i] - x2[i], 2);
    }
    return sqrt(norm);
  }

  std::vector< std::vector<int> > read_data(std::ifstream& input_file, const char delimiter) {
    /*

      Saves each row of data as a column in output_data
    
      The data is expected to be of the form:
      x1 y1 z1 ...
      x2 y2 z2 ...
      ...
      thus each "column" of output_data contains all
      input data of one variable

    */
    // determine number of columns by reading first line
    // count rows, not columns
    std::string line;
    int nrows = 0;
    while(std::getline(input_file, line)) {
      nrows++;
    }
    // move back to beginning of file
    // need to clear EOF flag first
    input_file.clear();
    input_file.seekg(0);
    std::string val;
    std::vector< std::vector<int> > output_data(nrows);
    int j = 0;
    while(std::getline(input_file, line)) {
      std::stringstream ss(line);
      while(std::getline(ss, val, delimiter)) {
	output_data[j].push_back(std::atoi(val.c_str()));
      }
      j++;
    }
    return output_data;
  }

} 
