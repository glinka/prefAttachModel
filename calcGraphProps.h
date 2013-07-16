#ifndef CALCGRAPHPROPS_H
#define CALCGRAPHPROPS_H

class calcGraphProps {
 private:
  calcGraphProps() {};
  ~calcGraphProps() {};
 public:
  static int *getDegrees(int **A, const int n);
  static int *getAdjEigVals(int **A, const int n);
  /**  static int **getAdjEigVects(int (&A)[n][m]);
  static int *getRandWalkLaplEigVals(int (&A)[n][m]);
  static int **getRandWalkLaplEigVects(int &A[n][m]);
  **/
};

#endif
