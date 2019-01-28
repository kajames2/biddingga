#include "numericaldists/line_ops.h"

#include <algorithm>
#include <cassert>
#include <vector>

#include "numericaldists/bounds.h"
#include "numericaldists/interval.h"
#include "numericaldists/line_segment.h"

namespace numericaldists {
std::vector<float> GetMesh(Interval interval, int n_points) {
  assert(n_points > 1);
  std::vector<float> vec(n_points);
  float seg_size = (interval.max - interval.min) / (n_points - 1);
  std::generate(vec.begin(), vec.end(), [n = interval.min, seg_size]() mutable {
    float res = n;
    n += seg_size;
    return res;
  });
  return vec;
}

std::vector<LineSegment> PointsToLines(const std::vector<float>& xs,
                                       const std::vector<float>& ys) {
  assert(xs.size() == ys.size());
  int n_lines = xs.size() - 1;
  std::vector<LineSegment> segments;
  segments.reserve(n_lines);
  for (int i = 0; i < n_lines; ++i) {
    segments.emplace_back(Bounds{xs[i], xs[i + 1]}, Bounds{ys[i], ys[i + 1]});
  }
  return segments;
}

}  // namespace numericaldists
