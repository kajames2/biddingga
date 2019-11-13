#ifndef BIDDINGGA_GRID_H_
#define BIDDINGGA_GRID_H_

#include <eigen3/Eigen/Core>

namespace numericaldists {

struct Grid {
  Eigen::ArrayXd xs;
  Eigen::ArrayXd ys;
  Eigen::ArrayXXd zs;
};

std::ostream& operator<<(std::ostream& os, const Grid& grid);
}  // namespace numericaldists

#endif  // BIDDINGGA_GRID_H_
