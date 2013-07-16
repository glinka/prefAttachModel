#ifndef CALCGRAPHPROPS_H
#define CALCGRAPHPROPS_H

class calcGraphProps {
 private:
  calcGraphProps() {};
  ~calcGraphProps() {};
 public:
  static int *getDegrees(int **A, const int n);
  static double *getAdjEigVals(int **A, const int n);
  static double **getAdjEigVects(int **A, const int n);
  static double *getLaplEigVals(int **A, const int n);
  static double **getLaplEigVects(int **A, const int n);
};

#endif
