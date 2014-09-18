#include <iostream>
#include <cstdlib>
#include <string>
#include <mpi.h>
#include "pamCPI.h"
#include "newton.h"
#include "gmres.h"
#include "newton_wrapper.h"

int main(int argc, char** argv) {
  // start MPI
  int mpierr = MPI_Init(NULL, NULL);
  if(mpierr != MPI_SUCCESS) {
    std::cout << "Error initializing MPI, terminating" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, mpierr);
  }
  int n = 100;
  int m = 10000;
  double kappa = 1.0;
  long int savetofile_interval = 1000;
  bool new_init = true;
  int proj_step = 50000;
  int collect_interval = 5000;
  int off_manifold_wait = 20000;
  int nmicrosteps = 100000;
  double h = 0.1;
  int nkrylov = 6;
  std::string init_type = "erdos";
  //parse command line args, could be moved to separate fn?
  for(int i = 1; i < argc; i++) {
    if(argv[i][0] == '-') {
      std::string current_label = argv[i];
      char *current_arg = argv[i+1];
      if(current_label == "-n" || current_label == "-nodes") {
	n = std::atoi(current_arg);
      }
      else if(current_label == "-m" || current_label == "-nEdges" || current_label == "-edges") {
	m = std::atoi(current_arg);
      }
      else if(current_label == "-k" || current_label == "-kappa") {
	kappa = std::atof(current_arg);
      }
      else if(current_label == "-save_interval" || current_label == "-savetofile_interval" || current_label == "-si") {
	savetofile_interval = (long int) (std::atof(current_arg) + 0.5); // parse_longint(current_arg);
      }
      else if(current_label == "-proj_step" || current_label == "-pStep" || current_label == "-projInterval") {
	proj_step = std::atoi(current_arg);
      }
      else if(current_label == "-collect_interval" || current_label == "-ci" || current_label == "-collectionInterval") {
	collect_interval = std::atoi(current_arg);
      }
      else if(current_label == "-off_manifold_wait" || current_label == "-omw") {
	off_manifold_wait = std::atoi(current_arg);
      }
      else if(current_label == "-nmicrosteps" || current_label == "-nms") {
	nmicrosteps = std::atoi(current_arg);
      }
      else if(current_label == "-init_type" || current_label == "-init") {
	init_type.assign(current_arg);
      }
      else if(current_label == "-withoutinit" || current_label == "-noinit") {
	new_init = false;
      }
      else if(current_label == "-h") {
	h = std::atof(current_arg);
      }
      else if(current_label == "-nkrylov") {
	nkrylov = std::atoi(current_arg);
      }
      else {
	std::cout << "Your entry of: " << current_label << " is an invalid argument" << std::endl;
      }
    }
  }

  if(!new_init) {
    // proj_step should be zero, as no
    // projection is taking place, but
    // if "proj_step < m" project_degrees()
    // will not run, and thus no data will be saved
    proj_step = 2*m;
    off_manifold_wait = 0;
  }
  int i;
  MPI_Comm_rank(MPI_COMM_WORLD, &i);
  double start_time;
  pamCPI model(n, m, kappa, proj_step, collect_interval, off_manifold_wait, nmicrosteps, savetofile_interval);
  if(init_type == "erdos") {
    model.initGraph();
  }
  else if(init_type == "complete") {
    model.init_complete_graph();
  }
  if(i == 0) {
    // time from root process
    start_time = time(NULL);
    // root process is responsible for running newton-gmres
    Newton newton(10, 1e-2, 100);
    Eigen::VectorXd x = newton.find_zero(newton_wrapper::F, n*Eigen::VectorXd::Ones(n), h, GMRES(1e-4, nkrylov), model);
    std::cout << "CELEBRAAAAAAAAAATE" << std::endl;
    // stop other processes from receiving input
    bool receive_data = false;
    MPI_Bcast(&receive_data, 1, MPI_INT, 0, MPI_COMM_WORLD);
    std::cout << x << std::endl;
  }
  else {
    // other processes listen for root's data which is used as
    // input to run_single_step(), and then send back results

    // WHY THE FUCK DOESN'T THIS VECTOR WORK
    // std::vector<int> deg_seq(n);
    int deg_seq[n];

    bool receive_data;
    MPI_Bcast(&receive_data, 1, MPI_INT, 0, MPI_COMM_WORLD);
    while(receive_data) {

      // MPI_Bcast(&deg_seq.front(), n, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(deg_seq, n, MPI_INT, 0, MPI_COMM_WORLD);

      // std::vector<int> new_deq_seq = model.run_single_step(deg_seq);
      std::vector<int> new_deq_seq = model.run_single_step(std::vector<int>(deg_seq, deg_seq + n));
   
      MPI_Bcast(&receive_data, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
  }
    
  if(i == 0) {
    double end_time = time(NULL);
    double elapsed_time = difftime(end_time, start_time);
    std::cout << "Wall time: " << elapsed_time << " s" << std::endl;
  }

  MPI_Finalize();
  // end MPI
  return 0;
}

