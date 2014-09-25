#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include "dmaps.h"
#include "gaussian_kernel.h"
#include "util_fns.h"
#include "embeddings.h"
#include "prefAttachModel.h"

int main(int argc, char** argv) {

  std::cout << "<---------------------------------------->" << std::endl;

  const int graph_size = 100;
  const int kappa = 1;
  const int ngraphs = std::atoi(argv[1]);
  const std::string init_type = "erdos";
  // run_interval calculated based on 5*n^3 total steps
  const int run_interval = 2*std::pow(graph_size, 3)/ngraphs;
  typedef std::vector< std::vector<int> > adj_mat;
  std::vector<adj_mat> pa_graphs(ngraphs);
  prefAttachModel model(graph_size, kappa);
  if(init_type == "erdos") {
    model.init_er_graph(std::pow(graph_size, 2));
  }
  else if(init_type == "complete") {
    model.init_complete_graph();
  }
  else {
    std::cout << "init type: " << init_type << " unrecognized, exiting." << std::endl;
    std::exit(1);
  }
  for(int i = 0; i < ngraphs; i++) {
    pa_graphs[i] = model.run_nsteps(run_interval);
  }

  // make evenly spaced number of spectral params (maybe should be logarithmically evenly spaced?)
  const int n_spec_params = 100;
  std::vector<double> spectral_params(n_spec_params);
  const double spec_param_max = 0.01;
  const double spec_param_min = 0.0001;
  const double dspec_param = (spec_param_max - spec_param_min)/(n_spec_params - 1);

  std::vector< std::vector<double> > graph_embeddings(ngraphs);
  for(int i = 0; i < n_spec_params; i++) {
    spectral_params[i] = spec_param_min + i*dspec_param;
  }
  for(int i = 0; i < ngraphs; i++) {
    graph_embeddings[i] = spectral_embedding(pa_graphs[i], spectral_params);
  }
  std::cout << "--> Graphs generated and embedded" << std::endl;
  std::cout << "--> Graph embeddings saved in: ./embedding_data" << std::endl;
  std::ofstream output_ges("./embedding_data/" + init_type + "_graph_embeddings.csv");
  utils::save_matrix(graph_embeddings, output_ges);
  output_ges.close();

  // DMAPS it
  const int k = std::atoi(argv[2]);
  std::vector<double> eigvals(k);
  std::vector< std::vector<double> > eigvects(k);
  std::vector< std::vector<double> > W(k);
  double median = utils::get_median(utils::get_squared_distances(graph_embeddings));
  std::cout << "--> Median squared distance: " << median << std::endl;
  Gaussian_Kernel gk(median);
  const int dmaps_success = dmaps::map(graph_embeddings, gk, eigvals, eigvects, W, k, 1e-12);
  if(dmaps_success != 1) {
    std::cout << "dmaps encountered an error" << std::endl;
  }
  std::cout << "--> DMAP computed" << std::endl;
  std::cout << "--> DMAP output saved in: ./embedding_data" << std::endl;

  std::ofstream output_eigvals("./embedding_data/dmaps_" + init_type + "_embedding_eigvals.csv");
  std::ofstream output_eigvects("./embedding_data/dmaps_" + init_type + "_embedding_eigvects.csv");

  utils::save_vector(eigvals, output_eigvals);
  utils::save_matrix(eigvects, output_eigvects);
  output_eigvals.close();
  output_eigvects.close();

  // DMAPS is hard, PCA it
  Eigen::MatrixXd X(ngraphs, n_spec_params);
  for(int i = 0; i < ngraphs; i++) {
    for(int j = 0; j < n_spec_params; j++) {
      X(i,j) = graph_embeddings[i][j];
    }
  }
  Eigen::MatrixXd V;
  Eigen::VectorXd S;
  const int pca_success = eigen_solver::pca(X, V, S, k);
  if(pca_success != 1) {
    std::cout << "pca encountered an error" << std::endl;
  }
  std::cout << "--> PCA completed" << std::endl;
  std::cout << "--> PCA output saved in: ./embedding_data" << std::endl;

  output_eigvals.open("./embedding_data/pca_" + init_type + "_embedding_eigvals.csv");
  output_eigvects.open("./embedding_data/pca_" + init_type + "_embedding_eigvects.csv");

  utils::save_vector(std::vector<double>(S.data(), S.data() + k), output_eigvals);
  for(int i = 0; i < k; i++) {
    eigvects[i] = std::vector<double>(V.data() + i*n_spec_params, V.data() + (i+1)*n_spec_params);
  }
  utils::save_matrix(eigvects, output_eigvects);
  output_eigvals.close();
  output_eigvects.close();

  std::cout << "<---------------------------------------->" << std::endl;
  
  return 0;
}
