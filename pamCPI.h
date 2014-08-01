#ifndef PAMCPI_H
#define PAMCPI_H
#include <string>
#include "prefAttachModel.h"

class pamCPI : public prefAttachModel {
 public:
  pamCPI(const int n, const int m, const double kappa, const int projStep, const int collectInterval, const int offManifoldWait, const int nMicroSteps, const int save_interval);
    ~pamCPI() {};
    void runCPI(const int nSteps, const std::string init_type, const std::string dir="./cpi_data_default/", const std::string run_id="");
 private:
    const int projStep, collectInterval, offManifoldWait, nMicroSteps, save_interval;
    void project(std::vector<std::vector<std::vector<double> > > &data, std::vector<double> &time, std::ofstream* projData, std::ofstream* eigVectData, std::ofstream* eigval_data);
};

#endif
