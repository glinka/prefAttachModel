#ifndef PREFATTACH_H
#define PREFATTACH_H
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
};

#endif
