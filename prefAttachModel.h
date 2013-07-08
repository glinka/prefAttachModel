#ifndef PREFATTACH_H
#define PREFATTACH_H
#include <random>

struct graphData {
    int *degSeq;
    int **A;
};

class prefAttachModel {
 private:
  const int n, m;
  const double kappa;
  double rnNormalization;
  int **A;
  int *degs;
  std::mt19937 *mt;
  double genURN();
  void initGraph();
  graphData step(bool saveFlag);
  int consistencyCheck();
  void saveData(graphData *data, int nSteps, int dataInterval);
 public:
  void run(int nSteps, int dataInterval);
  prefAttachModel(int n, int m, double kappa);
  ~prefAttachModel() {
      for(int i = 0; i < n; i++) {
	  delete[] A[i];
      }
      delete[] A;
      delete[] degs;
      delete mt;
  };
};

#endif

  

  
