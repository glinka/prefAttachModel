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
    void saveData(graphData *data, int nData, std::ofstream &fileHandle);
    void saveData(std::vector<std::vector<double> > &data, std::ofstream &fileHandle);
    void saveData(std::vector<double> &data, std::ofstream &fileHandle);
    void save_coeffs(const std::vector< std::vector< double > > &data, std::ofstream &fileHandle);
 public:
    void run(long int nSteps, int dataInterval);
    std::ofstream &createFile(std::string base, std::vector<double> &addtnlData, std::vector<std::string> &addtnlDataLabels);
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
