#ifndef _NUMERICALDISTS_PIECEWISE_LINEAR_
#define _NUMERICALDISTS_PIECEWISE_LINEAR_

#include <iostream>
#include <cassert>
#include <numeric>
#include <vector>

#include "numericaldists/line_ops.h"
#include "numericaldists/line_segment.h"
#include "numericaldists/uneven_piecewise_linear.h"

namespace numericaldists {

class PiecewiseLinear : public UnevenPiecewiseLinear {
 public:
  PiecewiseLinear() {}
  PiecewiseLinear(const std::vector<float>& ys, Interval interval)
      : segment_length_((interval.max - interval.min) / (ys.size() - 1)),
        min_(interval.min) {
    assert(ys.size() > 1);
    interval_ = interval;
    auto xs = GetMesh(interval, ys.size());
    segments_ = PointsToLines(xs, ys);
  }

  PiecewiseLinear(const std::vector<float>& xs, const std::vector<float>& ys)
      : UnevenPiecewiseLinear(xs, ys),
        segment_length_((xs.back() - xs.front()) / (xs.size() - 1)),
        min_(xs.front()) {}
  
  PiecewiseLinear(std::vector<LineSegment> segments)
      : UnevenPiecewiseLinear(segments) {
    min_ = segments_.front().GetXInterval().min;
    float max = segments_.back().GetXInterval().max;
    segment_length_ = (max - min_) / segments_.size();
  }
  
  LineSegment GetLine(float x) const override {    
    return segments_[(x - min_) / segment_length_];
  }

 private:
  float segment_length_;
  float min_;
};

}  // namespace numericaldists

#endif  // _NUMERICALDISTS_PIECEWISE_LINEAR_
