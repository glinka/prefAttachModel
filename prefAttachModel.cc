#include <cstdlib>
#include <iostream>
#include <random>
#include <chrono>
#include <sstream>
#include <iomanip>
#include "prefAttachModel.h"

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

prefAttachModel::prefAttachModel(const int n, const int m, const double kappa): m(m), kappa(kappa), n(n) {
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  mt = new mt19937(seed);
  rnNormalization = (double) (mt->max()+1);
};

double prefAttachModel::genURN() {
    return (*mt)()/(rnNormalization);
}

void prefAttachModel::init_complete_graph() {
  // init a complete graph with loops
  m = n*n/2;
  A = new int*[n];
  degs = new int[n];
  for(int i = 0; i < n; i++) {
    A[i] = new int[n];
  }
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      A[i][j] = 1;
    }
    degs[i] = n;
  }
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

/****************************** TODO ******************************
must recalc graph properties as in previous initGraph(), otherwise segfault
******************************* TODO ******************************/

void prefAttachModel::initGraph(int **newA) {
    int i, j;
    double max = 0, min = 100;
    m = 0;
    for(i = 0; i < n; i++) {
	degs[i] = 0;
	for(j = 0; j < n; j++) {
	    A[i][j] = newA[i][j];
	    degs[i] += A[i][j];
	    if(A[i][j] > max) {
		max = A[i][j];
	    }
	    if(A[i][j] < min) {
		min = A[i][j];
	    }
	}
	m += degs[i];
    }
    m /= 2;
    cout << max << "," << min << "\n";
}

graphData *prefAttachModel::step(bool saveFlag) {
  int i = 0, degCount = 0;
  int oldEdge = ((int) floor(2*m*genURN())) + 1;
  while(degCount < oldEdge) {
    degCount += degs[i++];
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
  while(sum <= p) {
    sum += (degs[i++]+kappa)/(2*m+n*kappa);
  }
  int uNew = --i;
  A[uOld][v]--;
  A[v][uOld]--;
  A[uNew][v]++;
  A[v][uNew]++;
  degs[uOld]--;
  degs[uNew]++;
  /**
  if(consistencyCheck() == 1) {
      cout << "step error" << endl;
  }
  **/
  if(saveFlag) {
      graphData *data = new graphData;
      data->degSeq = new int[n];
      data->A = new int*[n];
      for(i = 0; i < n; i++) {
	  data->A[i] = new int[n];
	  data->degSeq[i] = degs[i];
	  for(j = 0; j < n; j++) {
	      data->A[i][j] = A[i][j];
	  }
      }
      return data;
  }
}

ofstream* prefAttachModel::createFile(const string base, const string dir, vector<double> &addtnlData, vector<string> &addtnlDataLabels) {
  stringstream ss;
  ss << dir << base << "_" << n << "_" << m << "_" << kappa;
  for(unsigned int i = 0; i < addtnlData.size(); i++) {
      ss << "_" << addtnlData[i];
  }
  ss << ".csv";
  string fileName = ss.str();
  ofstream *file = new ofstream();
  file->open(fileName);
  *file << "n=" << n << ",";
  *file << "m=" << m << ",";
  *file << "kappa=" << kappa;
  for(unsigned int i = 0; i < addtnlData.size(); i++) {
      *file << "," << addtnlDataLabels[i] << "=" << addtnlData[i];
  }
  *file << "\n";
  return file;
}

void prefAttachModel::run(long int nSteps, long int dataInterval, string init_type) {
  graphData *data;
  //data will be appended to file every time SAVE_INTERVAL data pts are collected
  const int SAVE_INTERVAL = 1000;
  if(init_type == "erdos") {
    initGraph();
  }
  else if(init_type == "complete") {
    init_complete_graph();
  }
  // initGraph() is default
  else {
    initGraph();
  }
  vector<double> forFile;
  vector< vector<int> > degs_to_save(SAVE_INTERVAL);
  vector< long int > times_to_save(SAVE_INTERVAL);
  vector< double > densities_to_save(SAVE_INTERVAL);
  vector< double > selfloop_densities_to_save(SAVE_INTERVAL);
  forFile.push_back(dataInterval);
  vector<string> forFileStrs;
  forFileStrs.push_back("dataInterval");
  ofstream* paData = createFile("paData", "datadefault/", forFile, forFileStrs);
  ofstream* deg_data = createFile("deg_data", "datadefault/" , forFile, forFileStrs);
  ofstream* time_data = createFile("time_data", "datadefault/" , forFile, forFileStrs);
  ofstream* density_data = createFile("density_data", "datadefault/" , forFile, forFileStrs);
  ofstream* selfloop_density_data = createFile("selfloop_density_data", "datadefault/" , forFile, forFileStrs);
  //create filename and make header to csv
  int current_index = 0;
  for(long int i = 0; i < nSteps; i++) {
    if(i % dataInterval == 0) {
      graphData *d = step(true);
      degs_to_save[current_index] = vector<int>(d->degSeq, d->degSeq+n);
      // delete returned data
      for(int j = 0; j < n; j++) {
	delete[] d->A[j];
      }
      delete[] d->A;
      delete[] d->degSeq;
      delete d;
      times_to_save[current_index] = i;
      densities_to_save[current_index] = simplified_edge_density();
      selfloop_densities_to_save[current_index] = compute_selfloop_density();
      //save data every SAVE_INTERVAL times data is collected
      current_index++;
      if(current_index % SAVE_INTERVAL == 0) {
	// saveData(data, SAVE_INTERVAL, paData);
	  save_degrees(degs_to_save, *deg_data);
	  saveData<long int>(times_to_save, *time_data);
	  saveData< double >(densities_to_save, *density_data);
	  saveData< double >(selfloop_densities_to_save, *selfloop_density_data);
	  current_index = 0;
      }
    }
    else{
      step(false);
    }
  }
  paData->close();
  deg_data->close();
  time_data->close();
  delete paData;
  delete deg_data;
  delete time_data;
  delete density_data;
  delete selfloop_density_data;
  // delete &paData;
  //take out the garbage
}

template<typename T>
void prefAttachModel::saveData(vector < T > &data, ofstream &fileHandle) {
  int nData = data.size();
  for(int i = 0; i < nData; i++) {
    fileHandle << data[i] << endl;
  }
}
template void prefAttachModel::saveData<int>(vector < int > &data, ofstream &fileHandle);
template void prefAttachModel::saveData<double>(vector < double > &data, ofstream &fileHandle);

void prefAttachModel::save_degrees(const vector< vector<int> > &degs, ofstream &fileHandle) {
  vector<int>::const_iterator val;
  for(vector< vector<int> >::const_iterator v = degs.begin(); v != degs.end(); v++) {
    for(val = (*v).begin(); val != (*v).end() - 1; val++) {
      fileHandle << *val << ",";
    }
    fileHandle << *val << endl;
  }
}

void prefAttachModel::saveData(vector<vector<double> > &data, ofstream &fileHandle) {
    //assume data is nxn vector
    int tempDegs[n];
    int i, j;
    for(i = 0; i < n; i++) {
	tempDegs[i] = 0;
	for(j = 0; j < n; j++) {
	    tempDegs[i] += data[i][j];
	}
    }
    int sortedDegs[n][2];
    for(i = 0; i < n; i++) {
	sortedDegs[i][0] = tempDegs[i];
	sortedDegs[i][1] = i;
    }
    qsort(sortedDegs, n, 2*sizeof(int), compInt);
    int newIndex;
    int tempArray[n][n];
    for(i = 0; i < n; i++) {
	newIndex = sortedDegs[i][1];
	for(j = 0; j < n; j++) {
	    tempArray[i][j] = data[newIndex][j];
	}
    }
    for(j = 0; j < n; j++) {
	newIndex = sortedDegs[j][1];
	for(i = 0; i < n; i++) {
	    data[i][j] = tempArray[i][newIndex];
	}
    }
    for(i = 0; i < n; i++) {
      for(j = 0; j < n; j++) {
	fileHandle << data[i][j];
	if(j != n-1) {
	  fileHandle << ",";
	}
      }
      fileHandle << endl;
    }
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

void prefAttachModel::save_coeffs(const vector< vector< double > > &data, ofstream &fileHandle) {
  vector< double >::const_iterator val;
  for(vector< vector< double > >::const_iterator v = data.begin(); v != data.end(); v++) {
    for(val = (*v).begin(); val != (*v).end()-1; val++) {
      fileHandle << *val << ",";
    }
    fileHandle << *val << endl;
  }
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

double prefAttachModel::simplified_edge_density() {
  int zero_count = 0;
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      if(A[i][j] == 0) {
	zero_count++;
      }
    }
  }
  return (n*n - zero_count)/((double) n*n);
}

double prefAttachModel::compute_selfloop_density() {
  int selfloop_count = 0;
  for(int i = 0; i < n; i++) {
    if(A[i][i] > 0) {
      selfloop_count++;
    }
  }
  return ((double) selfloop_count)/n;
}
