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
 private:
    int m;
    const double kappa;
 protected:
    const int n;
    double rnNormalization;
    int **A;
    int *degs;
    std::mt19937 *mt;
    double genURN();
    void initGraph();
    void initGraph(int **newA);
    graphData step(bool saveFlag);
    int consistencyCheck();
    void save_degrees(const std::vector< std::vector<int> > &degs, std::ofstream &fileHandle);
    void saveData(graphData *data, int nData, std::ofstream &fileHandle);
    void saveData(std::vector<std::vector<double> > &data, std::ofstream &fileHandle);
    template <typename T>
      void saveData(std::vector< T > &data, std::ofstream &fileHandle);
    void save_coeffs(const std::vector< std::vector< double > > &data, std::ofstream &fileHandle);
 public:
    void run(long int nSteps, long int dataInterval);
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
