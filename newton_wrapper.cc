#include <algorithm>
#include <mpi.h>
#include "newton_wrapper.h"
#include "pamCPI.h"

// TESTING
#include <iostream>

namespace newton_wrapper {
  Eigen::VectorXd F(const Eigen::VectorXd& deg_seq, pamCPI model) {
  
    // ONLY RUN ON ROOT PROCESS

    // provides interface between pref attach model, Eigen and newton

    // returns x - \Phi(x) where 'x' is a degree sequence
    // and \Phi(x) is a single coarse step, taken with run_single_step()
    
    const int n = deg_seq.size();
    std::vector<int> deg_seq_copy(deg_seq.data(), deg_seq.data() + n);
    bool receive_data = true;
    MPI_Bcast(&receive_data, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&deg_seq_copy.front(), n, MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<int> projected_deg_seq_copy = model.run_single_step(deg_seq_copy);

    // sort both, and return difference of sorted sequences
    std::sort(deg_seq_copy.begin(), deg_seq_copy.end());
    std::sort(projected_deg_seq_copy.begin(), projected_deg_seq_copy.end());
    std::cout << "Degs span in: " << deg_seq_copy.front() << " " << deg_seq_copy.back() << std::endl;
    std::cout << "Degs span out: " << projected_deg_seq_copy.front() << " " << projected_deg_seq_copy.back() << std::endl;
    Eigen::VectorXd projected_deg_seq(n);
    Eigen::VectorXd sorted_deg_seq(n);
    for(int i = 0; i < n; i++) {
      projected_deg_seq(i) = projected_deg_seq_copy[i];
      sorted_deg_seq(i) = deg_seq_copy[i];
    }
    // std::cout << "1-norm of F(x): " << sorted_deg_seq.sum() - projected_deg_seq.sum() << std::endl;
    std::cout << "m: " << sorted_deg_seq.sum()/2 << std::endl;
    return sorted_deg_seq - projected_deg_seq;
  }
}
