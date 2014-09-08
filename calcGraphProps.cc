#include <algorithm>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "calcGraphProps.h"

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

double *calcGraphProps::getAdjEigVals(int **A, const int n) {
  int i, j;
  /**
     all of this only to realize that Map != Matrix, but it may be useful later
     to reference these very odd offsets that needed to be added to the stride
  int offset;
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
  Matrix<double, Dynamic, Dynamic>::Map<Matrix<double, Dynamic, Dynamic>, 0, Stride<Dynamic, Dynamic> > M(A[0], n, n, Stride<Dynamic, Dynamic>(1,n+offset));
  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) {
      if(A[i][j] != M(i,j)) {
	cout << "A != M" << endl;
      }
    }
  }
  **/
  Matrix<double, Dynamic, Dynamic> M(n,n);
  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) {
      M(i,j) = A[i][j];
    }
  }
  SelfAdjointEigenSolver<Matrix<double, Dynamic, Dynamic> > solver(n);
  solver.compute(M); //will also compute eigvects
  Matrix<double, Dynamic, 1> eigVals = solver.eigenvalues();
  double *eigValsOut = new double[n];
  for(i = 0; i < n; i++) {
    eigValsOut[i] = eigVals(i);
  }
  return eigValsOut;
}

double **calcGraphProps::getAdjEigVects(int **A, const int n) {
  int i, j;
  Matrix<double, Dynamic, Dynamic> M(n,n);
  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) {
      M(i,j) = A[i][j];
    }
  }
  SelfAdjointEigenSolver<Matrix<double, Dynamic, Dynamic> > solver(n);
  solver.compute(M); //will also compute eigvects
  /**
     instead of creating a Matrix on the heap, create multidimensional array,
     probably a bit more space efficient. otherwise, use below:
  Matrix<double, Dynamic, Dynamic> *eigVects = new Matrix<double, Dynamic, Dynamic>;
  *eigVects = solver.eigenvectors();
  **/
  Matrix<double, Dynamic, Dynamic> eigVects = solver.eigenvectors();
  /**
     visual confirmation of eigenvector calculation success:
  Matrix<double, Dynamic, Dynamic> eigVals = solver.eigenvalues();
  cout << M << endl;
  cout << eigVects*eigVals.asDiagonal()*eigVects.inverse() << endl;
  **/
  double **eigVectOut = new double*[n];
  for(i = 0; i < n; i++) {
    eigVectOut[i] = new double[n];
    for(j = 0; j < n; j++) {
      eigVectOut[i][j] = eigVects(i,j);
    }
  }
  return eigVectOut;
}

double *calcGraphProps::getLaplEigVals(int **A, const int n) {
    int i, j, deg;
  Matrix<double, Dynamic, Dynamic> adj(n,n);
  Matrix<double, Dynamic, Dynamic> degs = Matrix<double, Dynamic, Dynamic>::Zero(n,n);
  for(i = 0; i < n; i++) {
    deg = 0;
    for(j = 0; j < n; j++) {
      adj(i,j) = A[i][j];
      deg += A[i][j];
    }
    degs(i,i) = deg;
  }
  Matrix<double, Dynamic, Dynamic> laplacian = degs - adj;
  SelfAdjointEigenSolver<Matrix<double, Dynamic, Dynamic> > solver(n);
  solver.compute(laplacian); //will also compute eigvects
  Matrix<double, Dynamic, Dynamic> eigVals = solver.eigenvalues();
  double *eigValsOut = new double[n];
  for(i = 0; i < n; i++) {
    eigValsOut[i] = eigVals(i);
  }
  return eigValsOut;
}

double **calcGraphProps::getLaplEigVects(int **A, const int n) {
    int i, j, deg;
  Matrix<double, Dynamic, Dynamic> adj(n,n);
  Matrix<double, Dynamic, Dynamic> degs = Matrix<double, Dynamic, Dynamic>::Zero(n,n);
  for(i = 0; i < n; i++) {
    deg = 0;
    for(j = 0; j < n; j++) {
      adj(i,j) = A[i][j];
      deg += A[i][j];
    }
    degs(i,i) = deg;
  }
  Matrix<double, Dynamic, Dynamic> laplacian = degs - adj;
  SelfAdjointEigenSolver<Matrix<double, Dynamic, Dynamic> > solver(n);
  solver.compute(laplacian); //will also compute eigvects
  Matrix<double, Dynamic, Dynamic> eigVects = solver.eigenvectors();
  double **eigVectOut = new double*[n];
  for(i = 0; i < n; i++) {
    eigVectOut[i] = new double[n];
    for(j = 0; j < n; j++) {
      eigVectOut[i][j] = eigVects(i,j);
    }
  }
  return eigVectOut;
}

// returns sorted degrees
vector<int> calcGraphProps::get_sorted_degrees(int **A, const int n) {
  vector<int> degs(n, 0);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      degs[i] += A[i][j];
    }
  }
  std::sort(degs.begin(), degs.end());
  return degs;
}
