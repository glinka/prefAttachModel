#include <sstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <dmaps.h>
#include <util_fns.h>
#include <igraph.h>
#include "custom_util_fns.h"
#include "kernel_function.h"
#include "embeddings.h"
#include "prefAttachModel.h"
#include "calcGraphProps.h"

typedef std::vector< std::vector<double> > matrix;
typedef std::vector<double> vector;

int main(int argc, char** argv) {

  typedef std::vector< std::vector<int> > adj_mat;

  std::cout << "<---------------------------------------->" << std::endl;

  /* // set up spectral parameter range \lambda */
  /* const int n_spec_params = 10; // 100; */
  /* std::vector<double> spectral_params(n_spec_params); */
  /* const double spec_param_max = -1, spec_param_min = -1.5; // max, min log(lambda) */
  /* const double dspec_param = (spec_param_max - spec_param_min)/(n_spec_params - 1); */
  /* for(int i = 0; i < n_spec_params; i++) { */
  /*   spectral_params[i] = std::pow(10, spec_param_min + i*dspec_param); */
  /* } */
  std::vector< std::vector<double> > graph_embeddings;
  const int k = 6;
  // if argv[1] == 0, generate data
  // otherwise load pre-exisitng data
  if(std::atoi(argv[1]) == 0) {
      
    // define parameters concerning dataset generation
    const int graph_size = 50; // 100;
    const double log_kappa_min = 0, log_kappa_max = 2; // smallest and largest log(kappa) value, increment by 'dkappa'
    const int nkappas = 30;
    const double dlog_kappa = 1.0*(log_kappa_max - log_kappa_min)/nkappas;
    const int m_min = 500, m_max = 5000, nms = 30;
    const double dm = 1.0*(m_max - m_min)/nms; // fewest and most edges to add, with increment 'dm'
    const int nintervals = 1, interval_length = 2*graph_size*graph_size*graph_size/nintervals; // take a total of 2*n^3 steps, collecting a data point every 'interval_length' steps
    std::vector< std::vector<double> > rho_kappas(nintervals*nkappas*nms, std::vector<double>(2));
    std::vector< std::vector<int> > degree_seqs(nintervals*nkappas*nms);

    // set igraph constants
    const int nmotifs = 8;
    const int nmotifs_disconnected = 11;
    const double scaling = 100;
    graph_embeddings = std::vector< std::vector<double> >(nintervals*nkappas*nms, std::vector<double>(nmotifs));

    // run at different values of kappa and m, saving the embeddings
#pragma omp parallel for collapse(2)
    for(int i = 0; i < nkappas; i++) {
      for(int j = 0; j < nms; j++) {
	const double kappa = std::pow(10, log_kappa_min + i*dlog_kappa); // new kappa value, placed here in inner loop for OMP
	const int m = m_min + j*dm;
	prefAttachModel model(graph_size, m, kappa);
	model.init("erdos");
	for(int k = 0; k < nintervals; k++) {
	  std::vector< std::vector<int> > model_graph = model.run_nsteps(interval_length);
	  // graph_embeddings[(i*nms + j)*nintervals + k] = spectral_embedding(model_graph, spectral_params);
	  rho_kappas[(i*nms + j)*nintervals + k][0] = m;
	  rho_kappas[(i*nms + j)*nintervals + k][1] = kappa;
	  degree_seqs[(i*nms + j)*nintervals + k] = calcGraphProps::get_sorted_degrees(model_graph);

	  // count subgraphs with igraph
	  double edges[2*graph_size*graph_size];
	  int count = 0;
	  for(int ii = 0; ii < graph_size; ii++) {
	    for(int jj = ii+1; jj < graph_size; jj++) {
	      if(model_graph[ii][jj] > 0) {
		edges[count] = ii;
		edges[count + 1] = jj;
		count += 2;
	      }
	    }
	  }

	  igraph_vector_t igraph_edges;
	  igraph_vector_init_copy(&igraph_edges, edges, count);

	  igraph_t igraph_graph;
	  igraph_empty(&igraph_graph, graph_size, IGRAPH_UNDIRECTED);
	  igraph_add_edges(&igraph_graph, &igraph_edges, 0);

	  igraph_vector_t cut_probs; // wtf is this
	  igraph_vector_init(&cut_probs, nmotifs_disconnected); // wtf is this
	  igraph_vector_fill(&cut_probs, 0); // wtf, fill it w/ zeros
	  for(int ii = 3; ii < 5; ii++) {
	    igraph_vector_t motif_counts;
	    igraph_vector_init(&motif_counts, 0);
	    igraph_motifs_randesu(&igraph_graph, &motif_counts, ii, &cut_probs);
	    // only keep non-nan values which are as follows
	    // three vertex: [nan nan # #]
	    // four vertex: [nan nan nan nan # nan # # # # #]
	    if(ii == 3) {
	      graph_embeddings[(i*nms + j)*nintervals + k][0] = VECTOR(motif_counts)[2]/scaling;
	      graph_embeddings[(i*nms + j)*nintervals + k][1] = VECTOR(motif_counts)[3]/scaling;
	    }
	    else {
	      graph_embeddings[(i*nms + j)*nintervals + k][2] = VECTOR(motif_counts)[4]/scaling;
	      graph_embeddings[(i*nms + j)*nintervals + k][3] = VECTOR(motif_counts)[6]/scaling;
	      graph_embeddings[(i*nms + j)*nintervals + k][4] = VECTOR(motif_counts)[7]/scaling;
	      graph_embeddings[(i*nms + j)*nintervals + k][5] = VECTOR(motif_counts)[8]/scaling;
	      graph_embeddings[(i*nms + j)*nintervals + k][6] = VECTOR(motif_counts)[9]/scaling;
	      graph_embeddings[(i*nms + j)*nintervals + k][7] = VECTOR(motif_counts)[10]/scaling;
	    }
	  }

	}
      }
    }
    // scale to mean zero
    /* matrix scaled_embeddings = math_util_fns::scale_data(graph_embeddings); */
    // investigate effect of epsilon on kernel sums

    /* const int k = std::atoi(argv[1]); */
    /* std::stringstream ss(""); */
    /* ss << "k=" << k; */
    /* std::string header = ss.str(); */

    std::cout << "--> Graphs generated and embedded" << std::endl;
    std::cout << "--> Graph embeddings saved in: ./embedding_data/rho_kappa_graph_embeddings.csv" << std::endl;
    std::cout << "--> Graph parameters saved in: ./embedding_data/rho_kappa_params.csv" << std::endl;
    std::cout << "--> Graph deg. sequences saved in: ./embedding_data/rho_kappa_deg_seqs.csv" << std::endl;
    util_fns::save_matrix(graph_embeddings, "./embedding_data/rho_kappa_graph_embeddings.csv");
    util_fns::save_matrix(rho_kappas, "./embedding_data/rho_kappa_params.csv");
    util_fns::save_matrix(degree_seqs, "./embedding_data/rho_kappa_deg_seqs.csv");
  }
  else {
    graph_embeddings = util_fns::read_data("./embedding_data/rho_kappa_graph_embeddings.csv");
    std::cout << "--> Pre-existing graphs loaded from ./embedding_data/rho_kappa_graph_embeddings.csv" << std::endl;
  }
     

  /* // test test_kernels function over log-spaced epsilons */
  /* const int nkernels = 20; */
  /* const int lower_exp = -6, upper_exp = 4; */
  /* const double de = (upper_exp - lower_exp)/(nkernels - 1.0); // change in epsilon */
  /* vector epsilons(nkernels); */
  /* for (int i = 0; i < nkernels; i++) { */
  /*   epsilons[i] = std::pow(10, lower_exp + i*de); */
  /* } */
  /* std::vector<Kernel_Function> kernels; */
  /* for (int i = 0; i < nkernels; i++) { */
  /*   kernels.push_back(Kernel_Function(epsilons[i])); */
  /* } */
  /* std::vector<double> kernel_sums = dmaps::test_kernels(scaled_embeddings, kernels); */

  /* util_fns::save_vector(kernel_sums, "./kernel_sums.csv"); */
  /* util_fns::save_vector(epsilons, "./epsilons.csv"); */
  /* std::cout << "--> Kernel sums and epsilons saved in ./kernel_sums.csv and ./epsilons.csv" << std::endl; */

  double median = custom_utils::get_median(custom_utils::get_squared_distances(graph_embeddings));
  std::cout << "--> Median squared distance: " << median << std::endl;
  const double epsilon = std::atof(argv[2]); // 500;
  std::cout << "--> Actual epsilon value: " << epsilon << std::endl;
  
  // DMAPS it, epsilon \approx 0.01 appears reasonable
  std::vector<double> eigvals(k);
  std::vector< std::vector<double> > eigvects(k);
  Kernel_Function gk(epsilon); // Gaussian kernel
  const int dmaps_success = dmaps::map(graph_embeddings, gk, eigvals, eigvects, k, 1e-12);
  if(dmaps_success != 1) {
    std::cout << "dmaps encountered an error" << std::endl;
  }
  std::cout << "--> DMAP computed" << std::endl;
  std::cout << "--> DMAP output saved in: ./embedding_data" << std::endl;

  util_fns::save_vector(eigvals, "./embedding_data/dmaps_rho_kappa_embedding_eigvals.csv");
  util_fns::save_matrix(eigvects, "./embedding_data/dmaps_rho_kappa_embedding_eigvects.csv");
  // DMAPS is hard, PCA it
  // jk, pca is hard
  // Eigen::MatrixXd X(npts, n_spec_params);
  // for(int i = 0; i < npts; i++) {
  //   for(int j = 0; j < n_spec_params; j++) {
  //     X(i,j) = graph_embeddings[i][j];
  //   }
  // }
  // Eigen::MatrixXd V;
  // Eigen::VectorXd S;
  // const int pca_success = eigen_solver::pca(X, V, S, k);
  // if(pca_success != 1) {
  //   std::cout << "pca encountered an error" << std::endl;
  //   return 0;
  // }
  // std::cout << "--> PCA completed" << std::endl;
  // std::cout << "--> PCA output saved in: ./embedding_data" << std::endl;

  // output_eigvals.open("./embedding_data/pca_" + init_type + "_embedding_eigvals.csv");
  // output_eigvects.open("./embedding_data/pca_" + init_type + "_embedding_eigvects.csv");

  // utils::save_vector(std::vector<double>(S.data(), S.data() + k), output_eigvals);
  // for(int i = 0; i < k; i++) {
  //   eigvects[i] = std::vector<double>(V.data() + i*n_spec_params, V.data() + (i+1)*n_spec_params);
  // }
  // utils::save_matrix(eigvects, output_eigvects);
  // output_eigvals.close();
  // output_eigvects.close();

  std::cout << "<---------------------------------------->" << std::endl;
  
  return 0;
}
