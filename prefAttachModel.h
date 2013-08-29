#ifndef PREFATTACH_H
#define PREFATTACH_H
#include <random>
#include <fstream>

struct graphData {
    int *degSeq;
    int **A;
};

class prefAttachModel {
    //clean up what's protected, what's private
 private:
    const int n, m;
    const double kappa;
 protected:
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
 public:
    void run(long int nSteps, int dataInterval);
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
