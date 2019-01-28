#include "numericaldists/uneven_piecewise_linear.h"

#include <algorithm>
#include <boost/container_hash/hash.hpp>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <numeric>

#include "numericaldists/line_ops.h"
#include "numericaldists/bounds.h"
#include "numericaldists/interval.h"
#include "numericaldists/line_segment.h"

namespace numericaldists {

UnevenPiecewiseLinear::UnevenPiecewiseLinear(const std::vector<float>& xs,
                                             const std::vector<float>& ys)
    : segments_(PointsToLines(xs, ys)),
      interval_(xs.front(), xs.back()){}

UnevenPiecewiseLinear::UnevenPiecewiseLinear(
    std::vector<LineSegment> segments) {
  std::stable_sort(segments.begin(), segments.end(),
                   [](const LineSegment& l1, const LineSegment& l2) {
                     if (l1.GetXInterval().min != l2.GetXInterval().min) {
                       return l1.GetXInterval().min < l2.GetXInterval().min;
                     } else {
                       return l1.GetXInterval().max < l2.GetXInterval().max;
                     }
                   });
  interval_ = Interval{segments.front().GetXInterval().min,
                       segments.back().GetXInterval().max};
  segments_ = segments;
}

LineSegment UnevenPiecewiseLinear::GetLine(float x) const {
  return *std::find_if(
      segments_.begin(), segments_.end(),
      [x](const LineSegment& line) { return line.IsInInterval(x); });
}

float UnevenPiecewiseLinear::GetBid(float x) const {
  if (x <= interval_.min) {
    return segments_[0].GetYBounds().a;
  } else if (x >= interval_.max) {
    return segments_[segments_.size() - 1].GetYBounds().b;
  }
  auto line = GetLine(x);
  return line.GetY(x);
}

template <class T>
std::vector<int> GetOrderings(const std::vector<T>& vec) {
  std::vector<int> indices(vec.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(),
            [&vec](int i1, int i2) { return vec[i1] < vec[i2]; });
  return indices;
}

std::vector<float> UnevenPiecewiseLinear::GetBids(std::vector<float> xs) const {
  auto ordering = GetOrderings(xs);
  int line_ind = 0;
  std::vector<float> bids(xs.size());
  for (int ind : ordering) {
    float x = xs[ind];
    if (x <= interval_.min) {
      bids[ind] = segments_[0].GetYBounds().a;
    } else if (x >= interval_.max) {
      bids[ind] = segments_[segments_.size() - 1].GetYBounds().b;
    } else {
      while (!segments_[line_ind].IsInInterval(x)) {
        ++line_ind;
      }
      bids[ind] = segments_[line_ind].GetY(x);
    }
  }
  return bids;
}

bool operator==(const UnevenPiecewiseLinear& g1,
                const UnevenPiecewiseLinear& g2) {
  return g1.segments_ == g2.segments_;
}

std::size_t hash_value(const UnevenPiecewiseLinear& s) {
  boost::hash<std::vector<LineSegment>> hasher;
  return hasher(s.segments_);
}

std::ostream& operator<<(std::ostream& os, const UnevenPiecewiseLinear& pl) {
  const auto separator = "\n";
  const auto* sep = "";
  for (auto ls : pl.segments_) {
    os << sep << ls;
    sep = separator;
  }
  return os;
}

}  // namespace numericaldists
