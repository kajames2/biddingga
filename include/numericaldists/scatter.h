#ifndef BIDDINGGA_SCATTER_H_
#define BIDDINGGA_SCATTER_H_

#include <ostream>

#include <eigen3/Eigen/Core>

namespace numericaldists {

struct Scatter {
  Eigen::ArrayXd xs;
  Eigen::ArrayXd ys;
};

std::ostream& operator<<(std::ostream& os, const Scatter& points);

}  // namespace numericaldists

#endif  // BIDDINGGA_SCATTER_H_
