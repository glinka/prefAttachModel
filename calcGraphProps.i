%module calcGraphProps
%{
#include "calcGraphProps.h"
%}


%typemap(in) (int **A, int n) {
  if(PyTuple_Check($input)) {
    int n = PyTuple_Size($input);
    $2 = n;
    int i, j;
    $1 = new int*[n];
    for(i = 0; i < n; i++) {
      $1[i] = new int[n];
      for(j = 0; j < n; j++) {
	$1[i][j] = PyInt_AsLong(PyTuple_GetItem(PyTuple_GetItem($input, i), j));
      }
    }
  }
 }

%inline %{
  double getDblPtrVal(double **ptr, int i, int j) {
    return (double) ptr[i][j];
  }
  double getSnglPtrVal(double *ptr, int i) {
    return (double) ptr[i];
  }
%}



%typemap(freearg) (int **A, int n) {
  for(int i = 0; i < $2; i++) {
    delete[] $1[i];
  }
  delete[] $1;
 }

%include "calcGraphProps.h"
