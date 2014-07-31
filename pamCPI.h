#ifndef PAMCPI_H
#define PAMCPI_H
#include <vector>
#include "prefAttachModel.h"

class pamCPI : public prefAttachModel {
 public:
    pamCPI(const int n, const int m, const double kappa, const int projStep, const int collectInterval, const int offManifoldWait, const int nMicroSteps);
    ~pamCPI() {};
    void runCPI(const int nSteps);
 private:
    const int projStep, collectInterval, offManifoldWait, nMicroSteps;
    void project(std::vector<std::vector<std::vector<double> > > &data, std::vector<double> &time, std::ofstream* projData, std::ofstream* eigVectData);
};

#endif
