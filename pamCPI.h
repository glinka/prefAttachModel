#ifndef PAMCPI_H
#define PAMCPI_H
#include <string>
#include "prefAttachModel.h"

class pamCPI : public prefAttachModel {
 public:
  pamCPI(const int n, const int m, const double kappa, const int projStep, const int collectInterval, const int offManifoldWait, const int nMicroSteps, const int save_interval);
    ~pamCPI() {};
    void runCPI(const int nSteps, const std::string init_type, const std::string dir_base="./cpi_data_default/", const std::string run_id="", const bool new_init=true);
 private:
    const int projStep, collectInterval, offManifoldWait, nMicroSteps, save_interval;
    void project(std::vector<std::vector<std::vector<double> > > &data, std::vector<double> &time, std::ofstream* projData, std::ofstream* eigVectData, std::ofstream* eigval_data);
    std::vector<int> project_degs(const std::vector< std::vector<int> >& data, const std::vector<double>& times, const int proj_step, const std::string run_id, const std::string dir);
};

#endif
