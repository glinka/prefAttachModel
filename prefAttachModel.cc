#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <random>
#include <chrono>
#include "prefAttachModel.h"

using namespace std;

prefAttachModel::prefAttachModel(const int n, const int m, const double kappa): n(n), m(m), kappa(kappa) {
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  mt19937 mt(seed);
  rnNormalization = (double) mt.max();
  randomGenerator = mt;
};

double prefAttachModel::getURN() {
  return randomGenerator()/rnNormalization;
}

void prefAttachModel::initGraph() {
  int i, j;
  //init random number generator
  //init adjacency matrix and edge vector
  A = new int*[n];
  degs = new int[n]
  for(i = 0, i < n; i++) {
    A[i] = new int[n];
  }
  int nEdges = (n*(n-1))/2 + n;
  int *edges = new int[nEdges];
  int newEdge;
  //assign m edges uniformly
  for(i = 0; i < m; i++) {
    newEdge = (int) floor(nEdges*getURN());
    edges[newEdge] = edges[newEdge] + 1;
  }
  int index = 0;
  for(i = 0; i < n; i++) {
    for(j = i; j < n; j++) {
      A[i][j] = edges[index++];
      degs[i] = degs[i] + A[i][j];
      degs[j] = degs[j] + A[i][j];
    }
  }
}

graphData prefAttachModel::step() {
  int sum, uOld, uNew, v;
  int oldEdge = (int) floor(2*m*genURN());
  int i, j, degCount, edgeCount = 0;
  double p;
  while(degCount < oldEdge) {
    degCount += degs[i++];
  }
  oldEdge = (int) floor(degs[--i]*genURN());
  while(edgeCount < oldEdge) {
    edgeCount += A[i][j++];
  }
  if(genURN() > 0.5) {
    uOld = i;
    v = --j;
  }
  else {
    uOld = --j;
    v = i;
  }
  //above should be faster for large n, as
  //the commented block below must loop through
  //n^2 elements to find the old edge
  /**
    sum = 0;
    for(j = i; j < n; j++) {
      sum += A[i][j];
    }
    edgeCount += sum;
    i++;
  }
  edgeCount -= sum;
  i--;
  j = i;
  while(edgeCount < oldEdge) {
    edgeCount += A[i][j++];
  }
  j--;
  **/
  //find new edge based on linear preferential attachment
  p = genURN();
  sum = 0;
  i = 0;
  while(sum < p) {
    sum += (degs[i++]+kappa)/(2*m+n*k);
  }
  uNew = --i;
  A[uOld][v] = A[uOld][v] - 1;
  A[uNew][v] = A[uNew][v] + 1;
  graphData data;
  data.degSeq = degs;
  return data;
}

void prefAttachModel::run(int nSteps, int dataInterval) {
  graphData *data;
  data = new graphData[nSteps/dataInterval];
  initGraph();
  for(int i = 0; i < nSteps; i++) {
    if(i % dataInterval == 0) {
      data[i/dataInterval] = step();
    }
    else{
      step();
    }
  }
}
