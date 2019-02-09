#ifndef _NUMERICALDISTS_UNEVEN_PIECEWISE_LINEAR_
#define _NUMERICALDISTS_UNEVEN_PIECEWISE_LINEAR_

#include <vector>

#include "numericaldists/interval.h"
#include "numericaldists/bounds.h"

namespace numericaldists {

class UnevenPiecewiseLinear {
 public:
  UnevenPiecewiseLinear() {}
  UnevenPiecewiseLinear(std::vector<float> xs,
                        std::vector<float> ys);
  float GetBid(float x) const;
  std::vector<float> GetBids(std::vector<float> xs) const;
  float operator()(float x) const { return GetBid(x); }
  std::vector<float> operator()(std::vector<float> xs) const {
    return GetBids(xs);
  }

 private:
  std::vector<float> xs_;
  std::vector<float> ys_;
  Interval x_range_;
  Bounds y_range_;
};
}  // namespace numericaldists

#endif  // _NUMERICALDISTS_UNEVEN_PIECEWISE_LINEAR_
