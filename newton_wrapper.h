#ifndef NEWTON_WRAPPER_H_
#define NEWTON_WRAPPER_H_

#include <Eigen/Dense>
class pamCPI;

namespace newton_wrapper {
  Eigen::VectorXd F(const Eigen::VectorXd& deg_seq, pamCPI model);
}

#endif
