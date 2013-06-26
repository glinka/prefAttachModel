#include <stdlib.h>
#include "prefAttachModel.h"

using namespace std;

int main(int argc, char *argv[]) {
  int n = atoi(argv[1]);
  int m = atoi(argv[2]);
  double kappa = atof(argv[3]);
  int nSteps = atoi(argv[4]);
  int dataInterval = atoi(argv[5]);
  prefAttachModel m(n, m, kappa);
  m.run(nStpes, dataInterval);
}
