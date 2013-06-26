#ifndef PREFATTACH_H
#define PREFATTACH_H
#include <random>
class prefAttachModel {
 private:
  struct graphData {
    int *degSeq;
  };
  const int n, m;
  const double kappa;
  double rnNormalization;
  int **A;
  int *degs;
  std::mt19937 (*randomGenerator)();
  double getURN();
  void initGraph();
  graphData step();
 public:
  void run(int nSteps, int dataInterval);
  prefAttachModel(int n, int m, double kappa);
  ~prefAttachModel(){};
};

#endif

  

  
