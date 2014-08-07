#ifndef PREFATTACH_H
#define PREFATTACH_H
#include <chrono>
#include <random>
#include <fstream>
#include <string>
#include <vector>

struct graphData {
    int *degSeq;
    int **A;
};

class prefAttachModel {
    //clean up what's protected, what's private
 protected:
    int m;
    const double kappa;
    const int n;
    double rnNormalization;
    int **A;
    int *degs;
    std::mt19937 *mt;
    double genURN();
    void initGraph();
    void initGraph(int **newA);
    void init_complete_graph();
    void init_graph_loosehh(std::vector<int> degs);
    graphData *step(bool saveFlag);
    int consistencyCheck();
    double simplified_edge_density();
    double compute_selfloop_density();
    void save_degrees(const std::vector< std::vector<int> > &degs, std::ofstream &fileHandle);
    void saveData(graphData *data, int nData, std::ofstream &fileHandle);
    void saveData(std::vector<std::vector<double> > &data, std::ofstream &fileHandle);
    template <typename T>
      void saveData(const std::vector < T > &data, std::ofstream &fileHandle) {
      int nData = data.size();
      for(int i = 0; i < nData; i++) {
	fileHandle << data[i] << std::endl;
      }
    }
    void save_coeffs(const std::vector< std::vector< double > > &data, std::ofstream &fileHandle);
 public:
    void run(long int nSteps, long int dataInterval, std::string init_type);
    std::ofstream* createFile(const std::string base, const std::string dir, std::vector<double> &addtnlData, std::vector<std::string> &addtnlDataLabels);
    prefAttachModel(int n, int m, double kappa);
    ~prefAttachModel() {
	for(int i = 0; i < n; i++) {
	    delete[] A[i];
	}
	delete[] A;
	delete[] degs;
	delete mt;
    };

    // rule of three you dumbass
 prefAttachModel(const prefAttachModel& tocopy): m(tocopy.m), kappa(tocopy.kappa), n(tocopy.n) {
      A = new int*[n];
      degs = new int[n];
      for(int i = 0; i < n; i++) {
	degs[i] = tocopy.degs[i];
	A[i] = new int[n];
	for(int j = 0; j < n; j++) {
	  A[i][j] = tocopy.A[i][j];
	}
      }
      unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
      mt = new std::mt19937(seed);
      rnNormalization = (double) (mt->max()+1);
    }

    prefAttachModel& operator=(const prefAttachModel& rhs) {
      if(this == &rhs) {
	return *this;
      }
      // a wiser man would properly define some swap routine
      // or remove the 'const'-ness of the members altogether
      // to solve this problem of const-reassignment, but
      // I shall take the easier path and simply assume that
      // all such values are the same for both lhs and rhs
      else {
	for(int i = 0; i < n; i++) {
	  degs[i] = rhs.degs[i];
	  for(int j = 0; j < n; j++) {
	    A[i][j] = rhs.A[i][j];
	  }
	}
	m = rhs.m;
	return *this;
      }
    }

};

#endif
