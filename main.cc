#include <stdlib.h>
#include "prefAttachModel.h"

using namespace std;

int main(int argc, char *argv[]) {
  int n = 100;
  int m = 10000;
  double kappa = 0.5;
  int nSteps = 10000000;
  int dataInterval = 1000;
  int i;
  //parse command line args, could be moved to separate fn?
  for(i = 1; i < argc; i++) {
    if(argv[i][0] == '-') {
      string currentLabel = argv[i];
      char *currentArg = argv[i+1];
      if(currentLabel == "-n" || currentLabel == "-nodes") {
	n = atoi(currentArg);
      }
      else if(currentLabel == "-m" || currentLabel == "-nEdges" || currentLabel == "-edges") {
	m = atoi(currentArg);
      }
      else if(currentLabel == "-k" || currentLabel == "-kappa") {
	kappa = atof(currentArg);
      }
      else if(currentLabel == "-s" || currentLabel == "-nSteps") {
	nSteps = atoi(currentArg);
      }
      else if(currentLabel == "-di" || currentLabel == "-dataInterval") {
	dataInterval = atoi(currentArg);
      }
    }
  }
  prefAttachModel model(n, m, kappa);
  model.run(nSteps, dataInterval);
}
