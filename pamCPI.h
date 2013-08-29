#ifndef PAMCPI_H
#define PAMCPI_H
#include "prefAttachModel.h"

class pamCPI : public prefAttachModel {
 public:
    pamCPI(prefAttachModel *pam, const int projStep, const int collectInterval, const int offManifoldWait, const int nMicroSteps);
    ~pamCPI() {};
    void runCPI(const int nSteps);
 private:
    const prefAttachModel *pam;
    const int projStep, collectInterval, offManifoldWait, nMicroSteps;
};

#endif
