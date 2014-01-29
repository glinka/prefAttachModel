#include <iostream>
#include <stdlib.h>
#include "prefAttachModel.h"
#include <cstring>

using namespace std;
const int ASCII_CHAR_OFFSET = 48;

long int parse_longint(const char* number) {
  int ndigits = strlen(number);
  for(int i = 0; i < ndigits; i++) {
    long int base = 1;
    for(int j = 0; j < ndigits-i-1; j++) {
      base *= 10;
    }
    return base*(number[i] - ASCII_CHAR_OFFSET) + parse_longint(number + i + 1);
  }
}

int main(int argc, char *argv[]) {
  int n = 100;
  int m = 10000;
  double kappa = 0.5;
  long int nSteps = 10000000;
  long int dataInterval = 1000;
  int i;
  bool project = false;
  //CPI vars
  int projStep = 1000;
  int collectInterval = 1000;
  int offManifoldWait = 5000;
  int nMicroSteps = 20000;
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
	nSteps = (long int) (atof(currentArg) + 0.5);// parse_longint(currentArg);
      }
      else if(currentLabel == "-di" || currentLabel == "-dataInterval" || currentLabel == "-ci" || currentLabel == "-collection_interval") {
	dataInterval = (long int) (atof(currentArg) + 0.5);// parse_longint(currentArg);
      }
      else if(currentLabel == "-project") {
	project = true;
      }
      else if(currentLabel == "-projStep" || currentLabel == "-pStep" || currentLabel == "-projInterval") {
	project = true;
	projStep = atoi(currentArg);
      }
      else if(currentLabel == "-collectInterval" || currentLabel == "-ci" || currentLabel == "-collectionInterval") {
	project = true;
	collectInterval = atoi(currentArg);
      }
      else if(currentLabel == "-offManifoldWait" || currentLabel == "-omw") {
	project = true;
	offManifoldWait = atoi(currentArg);
      }
      else if(currentLabel == "-nMicroSteps" || currentLabel == "-nms") {
	project = true;
	nMicroSteps = atoi(currentArg);
      }
    }
  }
  if(project) {
    //    pamCPI model(n, m, kappa, projStep, collectInterval, offManifoldWait, nMicroSteps);
    //    model.runCPI(nSteps);
  }
  else {
     prefAttachModel model(n, m, kappa);
     model.run(nSteps, dataInterval);
  }
  return 0;
}
