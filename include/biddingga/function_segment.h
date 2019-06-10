#ifndef BIDDINGGA_FUNCTION_SEGMENT_H_
#define BIDDINGGA_FUNCTION_SEGMENT_H_

#include <eigen3/Eigen/Core>

namespace biddingga {

struct Scatter {
  Eigen::ArrayXd xs;
  Eigen::ArrayXd ys;
};

}  // namespace biddingga

#endif  // BIDDINGGA_FUNCTION_SEGMENT_H_
