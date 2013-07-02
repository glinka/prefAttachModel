#ifndef PREFATTACH_H
#define PREFATTACH_H
#include <random>

struct graphData {
  int *degSeq;
};

class prefAttachModel {
 private:
  const int n, m;
  const double kappa;
  double rnNormalization;
  int **A;
  int *degs;
  std::mt19937 mt;
  double genURN();
  void initGraph();
  graphData step();
 public:
  void run(int nSteps, int dataInterval);
  prefAttachModel(int n, int m, double kappa);
  ~prefAttachModel(){};
};

#endif

  

  
