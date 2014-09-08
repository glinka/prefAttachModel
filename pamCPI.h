#ifndef PAMCPI_H
#define PAMCPI_H
#include <string>
#include "prefAttachModel.h"

// TESTING
#include <iostream>

class pamCPI : public prefAttachModel {
 public:
  pamCPI(const int n, const int m, const double kappa, const int projStep, const int collectInterval, const int offManifoldWait, const int nMicroSteps, const int save_interval);
    ~pamCPI() {};
    void runCPI(const int nSteps, const std::string init_type, const std::string run_id="", const bool new_init=true);
    std::vector<int> run_single_step(const std::vector<int>& degree_seq);

 pamCPI(const pamCPI& tocopy): prefAttachModel(tocopy), projStep(tocopy.projStep), collectInterval(tocopy.collectInterval), offManifoldWait(tocopy.offManifoldWait), nMicroSteps(tocopy.nMicroSteps), save_interval(tocopy.save_interval) {
      /* std::cout << "here pamcopy" << std::endl; */
    }

    pamCPI& operator=(const pamCPI& rhs) {
      if(this == &rhs) {
	return *this;
      }
      else {
	for(int i = 0; i < n; i++) {
	  degs[i] = rhs.degs[i];
	  for(int j = 0; j < n; j++) {
	    A[i][j] = rhs.A[i][j];
	  }
	}
	m = rhs.m;
	/* std::cout << "here assgn" << std::endl; */
	return *this;
      }
    }

 private:
    const int projStep, collectInterval, offManifoldWait, nMicroSteps, save_interval;
    void project(std::vector<std::vector<std::vector<double> > > &data, std::vector<double> &time, std::ofstream* projData, std::ofstream* eigVectData, std::ofstream* eigval_data);
    std::vector<int> project_degs(const std::vector< std::vector<int> >& data, const std::vector<double>& times, int& proj_step, const std::string run_id, const std::string dir);
    
};

#endif
