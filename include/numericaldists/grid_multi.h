#ifndef BIDDINGGA_GRID_MULTI_H_
#define BIDDINGGA_GRID_MULTI_H_

#include <ostream>

#include <vector>
#include <eigen3/Eigen/Core>

namespace numericaldists {

struct GridMulti {
  Eigen::ArrayXd xs;
  Eigen::ArrayXd ys;
  std::vector<Eigen::ArrayXXd> z_sets;
};

std::ostream& operator<<(std::ostream& os, const GridMulti& grid);

}  // namespace numericaldists

#endif  // BIDDINGGA_GRID_MULTI_H_
