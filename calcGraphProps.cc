#include <iostream>
#include "calcGraphProps.h"
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

int *calcGraphProps::getDegrees(int **A, const int n) {
  int i, j;
  int *degs = new int[n];
  for(i = 0; i < n; i++) {
    degs[i] = 0;
    for(j = 0; j < n; j++) {
      degs[i] = degs[i] + A[i][j];
    }
  }
  return degs;
}

int *calcGraphProps::getAdjEigVals(int **A, const int n) {
  int i, j, offset;
  switch(n%4) {
  case 0:
    offset = 4;
    break;
  case 1:
    offset = 3;
    break;
  case 2:
    offset = 2;
    break;
  case 3:
    offset = 5;
    break;
  }
  Map<Matrix<int, Dynamic, Dynamic>, 0, Stride<Dynamic, Dynamic> > M(A[0], n, n, Stride<Dynamic, Dynamic>(1,n+offset));
  /**
     check offsets are correct
  **/
  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) {
      if(A[i][j] != M(i,j)) {
	cout << "A != M" << endl;
      }
    }
  }
}
  


