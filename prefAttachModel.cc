#include <cstdlib>
#include <iostream>
#include <random>
#include <chrono>
#include <sstream>
#include <iomanip>
#include "prefAttachModel.h"
#include "fitCurves.h"

using namespace std;

int compndArray(const void *c1v, const void *c2v) {
    double c1 = *(double *) c1v;
    double c2 = *(double *) c2v;
    /**
       if you wanted to sort based on a different column:
       const int (*c1)[n] = (int (*)[n]) c1v;
       as we're only interested in the first column, no need to cast
       to an array
    **/
    if(c1 > c2) {
	return 1;
    }
    else if(c2 > c1) {
	return -1;
    }
    else {
	return 0;
    }
}

int compInt(const void *p1, const void *p2) {
  //overflow be damned
  return( *(int *)p1 - *(int *)p2);
}

int compFlt(const void *p1, const void *p2) {
    double fp1 = *(double *) p1;
    double fp2 = *(double *) p2;
    if(fp1 > fp2) {
	return 1;
    }
    else if(fp2 > fp1) {
	return -1;
    }
    else {
	return 0;
    }
}

prefAttachModel::prefAttachModel(const int n, const int m, const double kappa): n(n), m(m), kappa(kappa) {
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  mt = new mt19937(seed);
  rnNormalization = (double) (mt->max()+1);
};

double prefAttachModel::genURN() {
    return (*mt)()/(rnNormalization);
}

void prefAttachModel::initGraph() {
  int i, j;
  //init random number generator
  //init adjacency matrix and edge vector
  A = new int*[n];
  degs = new int[n];
  for(i = 0; i < n; i++) {
    A[i] = new int[n];
  }
  int nEdges = (n*(n-1))/2 + n;
  int *edges = new int[nEdges];
  for(i = 0; i < nEdges; i++) {
      edges[i] = 0;
  }
  for(i = 0; i < n; i++) {
      degs[i] = 0;
  }
  int newEdge;
  //assign m edges uniformly
  for(i = 0; i < m; i++) {
    newEdge = (int) floor(nEdges*genURN());
    edges[newEdge] = edges[newEdge] + 1;
  }
  int index = 0;
  for(i = 0; i < n; i++) {
      A[i][i] = 2*edges[index++];
      degs[i] = degs[i] + A[i][i];
    for(j = i+1; j < n; j++) {
      A[i][j] = edges[index++];
      A[j][i] = A[i][j];
      degs[i] = degs[i] + A[i][j];
      degs[j] = degs[j] + A[j][i];
    }
  }
  /**
  if(consistencyCheck() == 1) {
      cout << "init error" << endl;
  }
  **/
  delete[] edges;
}

void prefAttachModel::initGraph(const int **newA) {
    int i, j;
    for(i = 0; i < n; i++) {
	for(j = 0; j < n; j++) {
	    A[i][j] = newA[i][j];
	}
    }
}

graphData prefAttachModel::step(bool saveFlag) {
  int i = 0, degCount = 0;
  int oldEdge = ((int) floor(2*m*genURN())) + 1;
  while(degCount < oldEdge) {
    degCount = degCount + degs[i++];
  }
  degCount -= degs[--i];
  int j = 0;
  while(degCount < oldEdge) {
      degCount += A[i][j++];
  }
  j--;
  int uOld, v;
  if(genURN() > 0.5) {
    uOld = i;
    v = j;
  }
  else {
    uOld = j;
    v = i;
  }
  //find new edge based on linear preferential attachment
  double p = genURN();
  double sum = 0;
  i = 0;
  while(sum < p) {
    sum += (degs[i++]+kappa)/(2*m+n*kappa);
  }
  int uNew = --i;
  A[uOld][v] = A[uOld][v] - 1;
  A[v][uOld] = A[v][uOld] - 1;
  A[uNew][v] = A[uNew][v] + 1;
  A[v][uNew] = A[v][uNew] + 1;
  degs[uOld] = degs[uOld] - 1;
  degs[uNew] = degs[uNew] + 1;
  /**
  if(consistencyCheck() == 1) {
      cout << "step error" << endl;
  }
  **/
  if(saveFlag) {
      graphData data;
      data.degSeq = new int[n];
      data.A = new int*[n];
      for(i = 0; i < n; i++) {
	  data.A[i] = new int[n];
	  data.degSeq[i] = degs[i];
	  for(j = 0; j < n; j++) {
	      data.A[i][j] = A[i][j];
	  }
      }
      return data;
  }
}

void prefAttachModel::run(long int nSteps, int dataInterval) {
  graphData *data;
  //data will be appended to file every time SAVE_INTERVAL data pts are collected
  const int SAVE_INTERVAL = 10;
  data = new graphData[SAVE_INTERVAL];
  initGraph();
  //create filename and make header to csv
  stringstream ss;
  ss << "paData_" << n << "_" << m << "_" << kappa << "_" << nSteps << "_" << dataInterval << ".csv";
  string fileName = ss.str();
  ofstream paData;
  paData.open(fileName);
  paData << "n=" << n << ",";
  paData << "m=" << m << ",";
  paData << "kappa=" << kappa << ",";
  paData << "nSteps=" << nSteps << ",";
  paData << "dataInterval=" << dataInterval << "\n";
  for(long int i = 0; i < nSteps; i++) {
    if((i+1) % dataInterval == 0) {
      int dataIndex = ((i+1) / dataInterval) % SAVE_INTERVAL;
      data[dataIndex] = step(true);
      //save data every SAVE_INTERVAL times data is collected
      if(dataIndex == 0) {
	  saveData(data, SAVE_INTERVAL, paData);
      }
    }
    else{
      step(false);
    }
  }
  paData.close();
  //take out the garbage
  for(int i = 0; i < SAVE_INTERVAL; i++) {
       for(int j = 0; j < n; j++) {
	  delete[] data[i].A[j];
      }
      delete[] data[i].A;
      delete[] data[i].degSeq;
  }
  delete[] data;
}

void prefAttachModel::saveData(graphData *data, int nData, ofstream &fileHandle) {
    int i, j, k;
    int sorted[nData][n+1][n+1];
    int toSort[n+1][n+1];
    //for each piece of data, fill
    //with n x n+1 array with the first column
    //composed of the degrees
    /**
       toSort[i] is organized as follows:

          0   deg(v1)  deg(v2)  ...  deg(vn)
       deg(v1)  A11      A12    ...    A1n
       deg(v2)  A21      A22    ...    A2n
          .      .        .      .      .
	  .      .        .      .      .
	  .      .        .      .      .
       deg(vn)  An1      An2    ...    Ann

       and is thus symmetric
    **/
    for(i = 0; i < nData; i++) {
	for(j = 0; j < n+1; j++) {
	    for(k = 0; k < n+1; k++) {
		if((j == 0) && (k == 0)) {
		    toSort[j][k] = 0;
		}
		else if(j == 0) {
		    toSort[j][k] = data[i].degSeq[k-1];
		}
		else if(k == 0) {
		    toSort[j][k] = data[i].degSeq[j-1];
		}
		else {
		    toSort[j][k] = data[i].A[j-1][k-1];
		}
	    }
	}
	//toSort[i] cannot be cast into a double pointer. because
	//it's non-contiguous memory. dick move, bjarne, dick move.
	int sortedDegs[n][2];
	for(j = 0; j < n; j++) {
	    sortedDegs[j][0] = toSort[0][j+1];
	    sortedDegs[j][1] = j;
	}
	qsort(sortedDegs, n, 2*sizeof(int), compInt);
	int newIndex;
	int tempArray[n+1][n];
	for(k = 0; k < n; k++) {
	  newIndex = sortedDegs[k][1] + 1;
	  for(j = 0; j < n+1; j++) {
	      tempArray[j][k] = toSort[j][newIndex];
	  }
	}
	for(j = 0; j < n; j++) {
	  newIndex = sortedDegs[j][1] + 1;
	  for(k = 0; k < n; k++) {
	      sorted[i][j+1][k+1] = tempArray[newIndex][k];
	  }
	}
	sorted[i][0][0] = 0;
	for(j = 0; j < n; j++) {
	    sorted[i][0][j+1] = sortedDegs[j][0];
	    sorted[i][j+1][0] = sortedDegs[j][0];
	}
    }
    /**
       output data into single csv file with following format:
       n, m, kappa, nSteps, dataInterval
       sortNow[0][j][k]
       sortNow[1][j][k]
       ...
    **/
    for(i = 0; i < nData; i++) {
	for(j = 0; j < n; j++) {
	    for(k = 0; k < n; k++) {
		fileHandle << (int) sorted[i][j+1][k+1];
		//don't append comma to last data pt
		if(k < n-1) {
		    fileHandle << ",";
		}
	    }
	    fileHandle << "\n";
	}
    }
    fileHandle.flush();
}



int prefAttachModel::consistencyCheck() {
    int sum, i , j;
    int edgeSum = 0;
    for(i = 0; i < n; i++) {
	sum = 0;
	for(j = 0; j < n; j++) {
	    sum += A[i][j];
	    if(A[i][j] != A[j][i]) {
		cout << "symmmetry broken" << endl;
		return 1;
	    }
	}
	edgeSum += sum;
	if(sum != degs[i]) {
	    cout << sum << " " << degs[i] << " " << i << endl;
	    return 1;
	}
	else if(sum < 0) {
	  cout << "negative degree at: " << i << " " << j << endl;
	}
    }
    if(edgeSum/2 != m) {
	cout << "nonconservative edge simulation" << endl;
	return 1;
    }
    return 0;
}
