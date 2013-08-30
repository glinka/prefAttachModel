#include <iostream>
#include <stdlib.h>
#include "pamCPI.h"

using namespace std;

int main(int argc, char *argv[]) {
  int n = 100;
  int m = 10000;
  double kappa = 0.5;
  long int nSteps = 10000000;
  int dataInterval = 1000;
  int i;
  //parse command line args, could be moved to separate fn?
  for(i = 1; i < argc; i++) {
    if(argv[i][0] == '-') {
      string currentLabel = argv[i];
      char *currentArg = argv[i+1];
      if(currentLabel == "-h" || currentLabel == "-help") {
	cout << "\nDefault parameters are:" << "\n";
	cout << "\tn : " << n << "\n";
	cout << "\tm : " << m << "\n";
	cout << "\tkappa: " << kappa << "\n";
	cout << "\tsteps: " << nSteps << "\n";
	cout << "\tcollection_interval: " << dataInterval << "\n";
	cout << "\nUse -variable_name to change values on the command line, e.g. \n./prefAttachModel -n 500 -m 250000 \nfor a 500 vertex, 250000 edge run with default kappa, steps and collection interval" << endl << endl;
	return 0;
      }
      else if(currentLabel == "-n" || currentLabel == "-nodes") {
	n = atoi(currentArg);
      }
      else if(currentLabel == "-m" || currentLabel == "-nEdges" || currentLabel == "-edges") {
	m = atoi(currentArg);
      }
      else if(currentLabel == "-k" || currentLabel == "-kappa") {
	kappa = atof(currentArg);
      }
      else if(currentLabel == "-s" || currentLabel == "-nSteps" || currentLabel == "-steps") {
	nSteps = atoi(currentArg);
      }
      else if(currentLabel == "-di" || currentLabel == "-dataInterval" || currentLabel == "-ci" || currentLabel == "-collection_interval") {
	dataInterval = atoi(currentArg);
      }
    }
  }
  int projStep = 10000;
  int collectInterval = 100;
  int offManifoldWait = 1000;
  int nMicroSteps = 2000;
  pamCPI model(n, m, kappa, projStep, collectInterval, offManifoldWait, nMicroSteps);
  model.runCPI(nSteps);
  /**
     prefAttachModel model(n, m, kappa);
     model.run(nSteps, dataInterval);
  **/
  return 0;
}
