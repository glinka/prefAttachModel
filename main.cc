#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <mpi.h>
#include "pamCPI.h"

using namespace std;
const int ASCII_CHAR_OFFSET = 48;

long int parse_longint(const char* number) {
  int ndigits = strlen(number);
  if(ndigits > 0) {
    long int base = 1;
    for(int j = 0; j < ndigits-1; j++) {
      base *= 10;
    }
    return base*(number[0] - ASCII_CHAR_OFFSET) + parse_longint(number + 1);
  }
  return 0;
}

string itos(const int i) {
  stringstream ss("");
  ss << i;
  return ss.str();
}

string create_dir(string dir_base) {
  stringstream ss;
  int folder_counter = 0;
  string dir;
    bool isdir = true;
    struct stat stat_dir;
    do {
      ss.str("");
      ss << dir_base << folder_counter << "/";
      folder_counter++;
      dir = ss.str();
      int check = stat(dir.c_str(), &stat_dir);
      if(check == -1) {
	mkdir(dir.c_str(), 0700);
	isdir = false;
      }
    } while (isdir);
    return dir;
}

int main(int argc, char *argv[]) {
  int n = 100;
  int m = 10000;
  double kappa = 0.5;
  long int nSteps = 10000000;
  long int savetofile_interval = 1000;
  int i;
  bool project = false;
  //CPI vars
  int projStep = 1000;
  int collectInterval = 1000;
  int offManifoldWait = 5000;
  int nMicroSteps = 20000;
  int nruns = 1;
  string init_type = "complete";
  int nthreads = 2;
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
	cout << "\tcollection_interval: " << savetofile_interval << "\n";
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
	nSteps = (long int) (atof(currentArg) + 0.5); // parse_longint(currentArg);
      }
      else if(currentLabel == "-save_interval" || currentLabel == "-savetofile_interval" || currentLabel == "-si") {
	savetofile_interval = (long int) (atof(currentArg) + 0.5); // parse_longint(currentArg);
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
      else if(currentLabel == "-init_type" || currentLabel == "-init") {
	init_type.assign(currentArg);
	cout << "initial graph is type: " << init_type << endl;
      }
      else if(currentLabel == "-nruns") {
	nruns = atoi(currentArg);
      }
      else if(currentLabel == "-nthreads") {
	nruns = atoi(currentArg);
      }
    }
  }
  if(project) {
    string dir = create_dir("./cpi_data");
    cout << "--> saving files into " << dir << endl;
    // #pragma omp parallel for num_threads(nthreads) schedule(dynamic)
    //     for(i = 0; i < nruns; i++) {
    // }

    // start MPI
    int mpierr = MPI_Init(NULL, NULL);
    if(mpierr != MPI_SUCCESS) {
      cout << "Error initializing MPI, terminating" << endl;
      MPI_Abort(MPI_COMM_WORLD, mpierr);
    }
    // end MPI

    pamCPI model(n, m, kappa, projStep, collectInterval, offManifoldWait, nMicroSteps, savetofile_interval);
    model.runCPI(nSteps, init_type, dir, itos(i));

    // start MPI
    MPI_Finalize();
    // end MPI

  }
  else {
    prefAttachModel model(n, m, kappa);
    for(i = 0; i < nruns; i++) {
      model.run(nSteps, savetofile_interval, init_type);
    }
  }
  return 0;
}
