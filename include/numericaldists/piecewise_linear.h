#ifndef _NUMERICALDISTS_PIECEWISE_LINEAR_
#define _NUMERICALDISTS_PIECEWISE_LINEAR_

#include <cassert>
#include <vector>

#include "numericaldists/interval.h"
#include "numericaldists/bounds.h"

namespace numericaldists {

class PiecewiseLinear {
 public:
  PiecewiseLinear() {}
  PiecewiseLinear(std::vector<float> ys, Interval interval);

  float GetBid(float x) const;
  std::vector<float> GetBids(std::vector<float> xs) const;

  float operator()(float x) const { return GetBid(x); }
  std::vector<float> operator()(std::vector<float> xs) const {
    return GetBids(xs);
  }

 private:
  std::vector<float> ys_;
  float inv_segment_length_;
  Interval x_range_;
  Bounds y_range_;
};

}  // namespace numericaldists

#endif  // _NUMERICALDISTS_PIECEWISE_LINEAR_
