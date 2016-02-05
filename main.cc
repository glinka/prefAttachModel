#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <mpi.h>
#include <util_fns.h>
#include "pamCPI.h"
#include "custom_util_fns.h"
#include "calcGraphProps.h"

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

// mpirun -n 100  ./pref_attach -n 100 -m 10000 -s 2000000 -project -offManifoldWait 10000 -nMicroSteps 60000 -save_interval 10000 -collectionInterval 1000 -kappa 1 -projStep 20000 // could do -projStep 30000
int main(int argc, char *argv[]) {
  int n = 100;
  int m = 10000;
  double kappa = 1;
  long int nSteps = 10000000;
  long int savetofile_interval = 1000;
  bool project = false;
  bool new_init = true;
  //CPI vars
  int projStep = 1000;
  int collectInterval = 1000;
  int offManifoldWait = 5000;
  int nMicroSteps = 20000;
  int nruns = 1;
  string init_type = "erdos";
  string input_filename = "./VSPECIAL_DATA/xs.csv";
  bool from_file = false;
  int nthreads = 2;
  //parse command line args, could be moved to separate fn?
  for(int i = 1; i < argc; i++) {
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
      }
      else if(currentLabel == "-withoutinit" || currentLabel == "-noinit") {
	new_init = false;
      }
      else if(currentLabel == "-input_filename") {
	input_filename.assign(currentArg);
	from_file = true;
      }
      else {
	cout << currentLabel << " not recognized as input argument" << endl;
      }
      // else if(currentLabel == "-nruns") {
      // 	nruns = atoi(currentArg);
      // }
      // else if(currentLabel == "-nthreads") {
      // 	nruns = atoi(currentArg);
      // }
    }
  }
  if(project) {
    if(!new_init) {
      // projStep should be zero, as no
      // projection is taking place, but
      // if "projStep < m" project_degrees()
      // will not run, and thus no data will be saved
      projStep = 2*m;
      offManifoldWait = 0;
    }
    // string dir = create_dir("./cpi_data");
    // cout << "--> saving files into " << dir << endl;
    // #pragma omp parallel for num_threads(nthreads) schedule(dynamic)
    //     for(i = 0; i < nruns; i++) {
    // }

    // start MPI //
    int mpierr = MPI_Init(NULL, NULL);
    if(mpierr != MPI_SUCCESS) {
      cout << "Error initializing MPI, terminating" << endl;
      MPI_Abort(MPI_COMM_WORLD, mpierr);
    }
    int i;
    MPI_Comm_rank(MPI_COMM_WORLD, &i);
    double start_time;
    if(i == 0) {
      // time from root process
      start_time = time(NULL);
    }

    pamCPI model(n, m, kappa, projStep, collectInterval, offManifoldWait, nMicroSteps, savetofile_interval);
    model.runCPI(nSteps, init_type, itos(i), new_init);

    if(i == 0) {
      double end_time = time(NULL);
      double elapsed_time = difftime(end_time, start_time);
      cout << "Wall time: " << elapsed_time << " s" << endl;
    }

    MPI_Finalize();
    // end MPI

  }
  else if(from_file) {
    // start MPI
    int mpierr = MPI_Init(NULL, NULL);
    if(mpierr != MPI_SUCCESS) {
      cout << "Error initializing MPI, terminating" << endl;
      MPI_Abort(MPI_COMM_WORLD, mpierr);
    }
    pamCPI model(n, m, kappa, projStep, collectInterval, offManifoldWait, nMicroSteps, savetofile_interval);
    model.run_fromfile(nSteps, input_filename);
    MPI_Finalize();
    // end MPI
  } 
  else {
    prefAttachModel model(n, kappa);
    if(init_type == "erdos") {
      model.init_er_graph(m);
    }
    else if(init_type == "complete") {
      model.init_complete_graph();
    }

    const int nintervals = 500;
    const int interval = nSteps/nintervals;
    ofstream degs_out("./pa_run_data/degs.out");
    ofstream times_out("./pa_run_data/times.out");
    vector< vector<int> > degs(nintervals);
    vector<int> times(nintervals);

    for(int i = 0; i < nintervals; i++) {
      /* degs[i] = calcGraphProps::get_sorted_degrees(model.run_nsteps(interval)); */
      times[i] = (i+1)*interval;
    }

    util_fns::save_matrix(degs, degs_out);
    util_fns::save_vector(times, times_out);
  }
  return 0;
}
