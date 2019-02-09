#include "numericaldists/uneven_piecewise_linear.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

namespace numericaldists {

UnevenPiecewiseLinear::UnevenPiecewiseLinear(std::vector<float> xs,
                                             std::vector<float> ys)
    : xs_(std::move(xs)),
      ys_(std::move(ys)),
      x_range_(xs_.front(), xs_.back()),
      y_range_(ys_.front(), ys_.back()) {
  assert(xs_.size() > 1);
  assert(ys_.size() == xs_.size());
  assert(std::is_sorted(xs_.begin(), xs_.end()));
}

float UnevenPiecewiseLinear::GetBid(float x) const {
  if (x <= x_range_.min) {
    return y_range_.a;
  }
  if (x >= x_range_.max) {
    return y_range_.b;
  }
  int i = 1;
  while (x > xs_[i]) {
    ++i;
    assert(i <= xs_.size());
  }
  float alpha = (x - xs_[i - 1]) / (xs_[i] - xs_[i - 1]);
  return (1 - alpha) * ys_[i - 1] + alpha * ys_[i];
}

std::vector<float> UnevenPiecewiseLinear::GetBids(std::vector<float> xs) const {
  assert(std::is_sorted(xs.begin(), xs.end()));
  std::vector<float> res(xs.size());
  int start_i = 0;
  int end_i = xs.size() - 1;
  while (start_i <= end_i && xs[start_i] <= x_range_.min) {
    res[start_i] = y_range_.a;
    ++start_i;
  }
  while (end_i >= start_i && xs[end_i] >= x_range_.max) {
    res[end_i] = y_range_.b;
    --end_i;
  }
  int n = 1;
  for (int i = start_i; i <= end_i; ++i) {
    while (xs[i] > xs_[n]) {
      ++n;
    }
    float alpha = (xs[i] - xs_[n - 1]) / (xs_[n] - xs_[n - 1]);
    res[i] =  (1 - alpha) * ys_[n - 1] + alpha * ys_[n];
  }
  return res;
}

}  // namespace numericaldists
