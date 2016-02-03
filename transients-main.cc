#include <vector>
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
#include "custom_util_fns.h"
#include "calcGraphProps.h"
#include "prefAttachModel.h"

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
  bool from_file = false;
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
  prefAttachModel model(n, m, kappa);
  vector<int> tri_init_deg_seq = model.init_triangle_graph();
  
  const int nintervals = 500;
  const int interval = nSteps/nintervals;

  vector<int> tri_tris(nintervals+1);
  tri_tris[0] = model.count_triangles();

  for(int i = 0; i < nintervals; i++) {
    model.run_nsteps(interval);
    tri_tris[i+1] = model.count_triangles();
  }

  ofstream tri_tris_out("./paper-data/other-runs/tri_tris.csv");
  util_fns::save_vector(tri_tris, tri_tris_out);

  // now restart with rando
  model.init_rando_graph(tri_init_deg_seq);

  vector<int> rando_tris(nintervals+1);
  rando_tris[0] = model.count_triangles();
  
  for(int i = 0; i < nintervals; i++) {
    model.run_nsteps(interval);
    rando_tris[i+1] = model.count_triangles();
  }

  ofstream rando_tris_out("./paper-data/other-runs/rando_tris.csv");
  util_fns::save_vector(rando_tris, rando_tris_out);


  // now restart with loosehh
  model.init_graph_loosehh(tri_init_deg_seq);

  vector<int> hh_tris(nintervals+1);
  hh_tris[0] = model.count_triangles();
  vector<int> times(nintervals+1);
  times[0] = 0;
  
  for(int i = 0; i < nintervals; i++) {
    model.run_nsteps(interval);
    hh_tris[i+1] = model.count_triangles();
    times[i+1] = (i+1)*interval;
  }

  ofstream hh_tris_out("./paper-data/other-runs/hh_tris.csv");
  util_fns::save_vector(hh_tris, hh_tris_out);
  ofstream times_out("./paper-data/other-runs/tris_times.csv");
  util_fns::save_vector(times, times_out);

  return 0;
}
