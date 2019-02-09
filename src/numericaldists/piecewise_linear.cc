#include "numericaldists/piecewise_linear.h"

#include <algorithm>
#include <cassert>

namespace numericaldists {

PiecewiseLinear::PiecewiseLinear(std::vector<float> ys, Interval interval)
    : ys_(std::move(ys)),
      inv_segment_length_((ys_.size() - 1) / (interval.max - interval.min)),
      x_range_(interval),
      y_range_(ys_.front(), ys_.back()) {
  assert(ys_.size() > 1);
}

float PiecewiseLinear::GetBid(float x) const {
  float x_adj = x - x_range_.min;
  if (x_adj < 0) {
    return y_range_.a;
  }
  if (x >= x_range_.max) {
    return y_range_.b;
  }
  float n_float = x_adj * inv_segment_length_;
  int n = static_cast<int>(n_float);
  float alpha = n_float - n;
  return (1 - alpha) * ys_[n] + alpha * ys_[n + 1];
}

std::vector<float> PiecewiseLinear::GetBids(std::vector<float> xs) const {
  assert(std::is_sorted(xs.begin(), xs.end()));
  int start_i = 0;
  int end_i = xs.size() - 1;
  std::vector<float> x_adjs(xs.size());
  std::transform(xs.begin(), xs.end(), x_adjs.begin(), [this](float x) {return x - x_range_.min;});
  std::vector<float> res(xs.size());
  while (start_i <= end_i && x_adjs[start_i] <= 0) {
    res[start_i] = y_range_.a;
    ++start_i;
  }
  while (end_i >= start_i && xs[end_i] >= x_range_.max) {
    res[end_i] = y_range_.b;
    --end_i;
  }

  for (int i = start_i; i <= end_i; ++i) {
    float n_float = x_adjs[i] * inv_segment_length_;
    int n = static_cast<int>(n_float);
    float alpha = n_float - n;
    res[i] = (1 - alpha) * ys_[n] + alpha * ys_[n + 1];
  }
  return res;
}

}  // namespace numericaldists
