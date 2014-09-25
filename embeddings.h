#ifndef EMBEDDINGS_FN_H
#define EMBEDDINGS_FN_H
#include <algorithm>
#include <vector>
/* #include <cmath> */
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
/* #include <Eigen/Geometry> */

std::vector<double> spectral_embedding(const Eigen::MatrixXd& G, const std::vector<double>& params, const std::vector<double>& start_dist, const std::vector<double>& stop_dist) {

  const int n = G.cols();

  // follow notation from paper
  Eigen::VectorXd q(n);
  Eigen::VectorXd p(n);
  for(int i = 0; i < n; i++) {
    q(i) = start_dist[i];
    p(i) = stop_dist[i];
  }

  /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
  /* !!!!!!!!!!!!!!!!!!!! Would (almost certainly) be !!!!!!!!!!!!!!!!!!!! */
  /* !!!!!!!!!!!!!!!!!!!! to search for eigenvalues that are the same !!!!!!!!!!!!!!!!!!!! */
  /* !!!!!!!!!!!!!!!!!!!! and Gram-Schmidt them to obtain P_inv !!!!!!!!!!!!!!!!!!!! */
  /* !!!!!!!!!!!!!!!!!!!! instead of naively, directly computing it !!!!!!!!!!!!!!!!!!!! */
  /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

  // assume undirected graph
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(G);
  Eigen::VectorXd D = es.eigenvalues();
  Eigen::MatrixXd P = es.eigenvectors();
  Eigen::VectorXd l_T = q.transpose()*P;
  Eigen::VectorXd r = P.inverse()*p;
  
  const int nparams = params.size();
  std::vector<double> S(nparams);
  for(int i = 0; i < nparams; i++) {
    S[i] = l_T.dot(Eigen::VectorXd(((D*params[i]).array().exp()*r.array())));
  }
  return S;
}

std::vector<double> spectral_embedding(const Eigen::MatrixXd& G, const std::vector<double>& params) {
  // cheap way to allow default arguments with const parameters
  // default to uniform starting and
  // stopping distributions
  const int n = G.cols();
  std::vector<double> start_dist(n, 1.0/n);
  std::vector<double> stop_dist(n, 1.0/n);
  return spectral_embedding(G, params, start_dist, stop_dist);
}  
  
std::vector<double> spectral_embedding(const std::vector< std::vector<int> >& G, const std::vector<double>& params) {
  const int n = G.size();
  Eigen::MatrixXd G_eigen(n, n);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      G_eigen(i,j) = G[i][j];
    }
  }
  return spectral_embedding(G_eigen, params);
}

#endif
