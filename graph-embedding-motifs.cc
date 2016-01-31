#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <util_fns.h>
#include <dmaps.h>
#include <kernel_function.h>
#include <igraph.h>
#include "custom_util_fns.h"
#include "embeddings.h"
#include "prefAttachModel.h"
#include "calcGraphProps.h"

int main(int argc, char** argv) {

  typedef std::vector< std::vector<int> > adj_mat;

  std::cout << "<---------------------------------------->" << std::endl;

  const int graph_size = 200;
  const int kappa = 1;
  // const int m = graph_size*(graph_size + 1)/2;
  // prefAttachModel model(graph_size, m, kappa);

  const int nintervals = std::atoi(argv[1]);
  const int run_interval = std::atoi(argv[2]);
  // run_interval calculated based on 5*n^3 total steps
  // const int run_interval = 6*std::pow(graph_size, 3)/ngraphs_per_type;
  // const int run_interval = 6*std::pow(graph_size, 2)/ngraphs_per_type;

  std::vector<std::string> init_types;
  // init_types.push_back("erdos");
  // init_types.push_back("erdos");
  // init_types.push_back("erdos");
  // init_types.push_back("erdos");
  // init_types.push_back("complete");
  // init_types.push_back("complete");
  // init_types.push_back("complete");
  // init_types.push_back("complete");

  // // new best 01/2016
  init_types.push_back("erdos");
  init_types.push_back("lopsided");
  // run with ./graph_embedding 800 10000 5 0.5 (last number can range in (0.1, 1)
  
  // // workable parameters:
  // // compile with
  // init_types.push_back("erdos");
  // init_types.push_back("complete");
  // // run with
  // ./graph_embedding 400 400 5
  // init_types.push_back("lopsided");
  // // better parameters:
  // // compile with
  /* init_types.push_back("erdos"); */
  /* init_types.push_back("complete"); */
  // // run with
  // ./graph_embedding 2000 40 5
  // // also good parameters:
  // // compile with
  // init_types.push_back("erdos");
  // init_types.push_back("erdos");
  // init_types.push_back("erdos");
  // init_types.push_back("erdos");
  // // run with
  // ./graph_embedding 100 10 5

  // init_types.push_back("lopsided");
  const int ntypes = init_types.size();
  const int npts = nintervals*ntypes;
  std::vector<adj_mat> pa_graphs(npts);
  std::vector< std::vector<int> > pa_degs(npts);
  std::vector<int> pa_times(npts);
  for(int i = 0; i < ntypes; i++) {
    // initialize with some fraction of nChoose2 edges
    // int m = (1.0*(i+1)/ntypes)*graph_size*(graph_size + 1)/2;
    // initialize with 50% of nChoose2 edges
    int m = 0.5*graph_size*(graph_size + 1)/2;
    prefAttachModel model(graph_size, m, kappa);
    model.init(init_types[i]);
    for(int j = 0; j < nintervals; j++) {
      pa_graphs[i*nintervals+j] = model.run_nsteps(run_interval);
      pa_degs[i*nintervals+j] = calcGraphProps::get_sorted_degrees(pa_graphs[i*nintervals+j]);
      pa_times[i*nintervals+j] = (j+1)*run_interval;
    }
  }

  const std::string init_type = "many";
  // save degree sequences and times
  std::ofstream output_ds("./embedding_data/" + init_type + "_degs.csv");
  std::ofstream output_times("./embedding_data/" + init_type + "_times.csv");
  output_ds << "ntypes=" << ntypes << std::endl;
  output_times << "ntypes=" << ntypes << std::endl;
  util_fns::save_matrix(pa_degs, output_ds);
  util_fns::save_vector(pa_times, output_times);
  output_ds.close();
  std::cout << "--> Graphs generated" << std::endl;

  // translate to igraph-useable format, count subraphs up to 4 vertices
  const int nmotifs = 9; // should be 9 diff simple motifs of <=4 vertices
  const int nthreads = 4;
  std::vector< std::vector<double> > graph_embeddings(npts, std::vector<double>(nmotifs));
  #pragma omp parallel for num_threads(nthreads) schedule(dynamic)
  for(int k = 0; k < npts; k++) {
    double edges[2*graph_size*graph_size];
    int count = 0;
    for(int i = 0; i < graph_size; i++) {
      for(int j = i+1; j < graph_size; j++) {
	if(pa_graphs[k][i][j] >= 1) {
	  edges[count] = i;
	  edges[count + 1] = j;
	  count += 2;
	}
      }
    }

    igraph_vector_t igraph_edges;
    igraph_vector_init_copy(&igraph_edges, edges, count);

    igraph_t igraph;
    igraph_empty(&igraph, graph_size, false);
    igraph_add_edges(&igraph, &igraph_edges, 0);

    igraph_vector_t cut_probs; // wtf is this
    igraph_vector_init(&cut_probs, 9); // wtf is this
    igraph_vector_fill(&cut_probs, 0); // wtf, fill it w/ zeros
    igraph_vector_t motif_counts;
    igraph_vector_init(&motif_counts, 0);
    igraph_motifs_randesu(&igraph, &motif_counts, 4, &cut_probs);
    // get back to the nicer vector
    for(int i = 0; i < nmotifs; i++) {
      graph_embeddings[k][i] = (int) VECTOR(motif_counts)[i];
    }
    std::cout << "--> Finished embedding " << k << " of " << npts << " graphs..." << std::endl;
  }
  
  std::cout << "--> Graphs embedded" << std::endl;
  std::cout << "--> Graph embeddings saved in: ./embedding_data" << std::endl;
  std::ofstream output_ges("./embedding_data/" + init_type + "_graph_embeddings.csv");
  util_fns::save_matrix(graph_embeddings, output_ges);
  output_ges.close();

  /* graph_embeddings = util_fns::read_data("./embedding_data/many_graph_embeddings.csv"); */
  /* std::cout << "--> Graphs loaded from ./embedding_data/many_graph_embeddings.csv" << std::endl; */

  // DMAPS it
  const int k = std::atoi(argv[3]);
  std::vector<double> eigvals(k);
  std::vector< std::vector<double> > eigvects(k);
  std::vector< std::vector<double> > W(k);
  double median = custom_utils::get_median(custom_utils::get_squared_distances(graph_embeddings));
  std::cout << "--> Median squared distance: " << median << std::endl;
  median = std::atof(argv[4]); // 0.8 6.0;
  std::cout << "--> Actual epsilon value in use: " << median << std::endl;
  Kernel_Function gk(median);
  const int dmaps_success = dmaps::map(graph_embeddings, gk, eigvals, eigvects, W, k, 1e-12);
  if(dmaps_success != 1) {
    std::cout << "dmaps encountered an error" << std::endl;
  }
  std::cout << "--> DMAP computed" << std::endl;
  std::cout << "--> DMAP output saved in: ./embedding_data" << std::endl;

  std::ofstream output_eigvals("./embedding_data/dmaps_" + init_type + "_embedding_eigvals.csv");
  std::ofstream output_eigvects("./embedding_data/dmaps_" + init_type + "_embedding_eigvects.csv");

  output_eigvals << "ntypes=" << ntypes << std::endl;
  output_eigvects << "ntypes=" << ntypes << std::endl;
  util_fns::save_vector(eigvals, output_eigvals);
  util_fns::save_matrix(eigvects, output_eigvects);
  output_eigvals.close();
  output_eigvects.close();

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

  // util_fns::save_vector(std::vector<double>(S.data(), S.data() + k), output_eigvals);
  // for(int i = 0; i < k; i++) {
  //   eigvects[i] = std::vector<double>(V.data() + i*n_spec_params, V.data() + (i+1)*n_spec_params);
  // }
  // util_fns::save_matrix(eigvects, output_eigvects);
  // output_eigvals.close();
  // output_eigvects.close();

  std::cout << "<---------------------------------------->" << std::endl;
  
  return 0;
}
