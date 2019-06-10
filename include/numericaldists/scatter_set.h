#ifndef NUMERICALDISTS_SCATTER_SET_H_
#define NUMERICALDISTS_SCATTER_SET_H_

#include <eigen3/Eigen/Core>

namespace numericaldists {

struct ScatterSet {
  Eigen::ArrayXd xs;
  Eigen::ArrayXXd ys;
};

}  // namespace numericaldists

#endif  // NUMERICALDISTS_SCATTER_SET_H_
