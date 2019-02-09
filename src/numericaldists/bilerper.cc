#include "numericaldists/bilerper.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

#include "numericaldists/piecewise_linear.h"

namespace numericaldists {

Bilerper::Bilerper(Interval x_int, Interval y_int,
                   std::vector<std::vector<float>> zs)
    : y_spacing_(GetSpan(y_int) / (zs.size() - 1)), y_int_(y_int) {
  for (auto slice = zs.rbegin(); slice != zs.rend(); ++slice) {
    slices_.push_back(PiecewiseLinear(*slice, x_int));
  }
}

float Bilerper::GetBid(float x, float y) const {
  if (y < y_int_.min) {
    return slices_[0](x);
  } else if (y >= y_int_.max) {
    return slices_[slices_.size() - 1](x);
  } else {
    float n_float = (y - y_int_.min) / y_spacing_;
    int n = static_cast<int>(n_float);
    float alpha = n_float - n;
    return (1 - alpha) * slices_[n](x) + alpha * slices_[n + 1](x);
  }
}

std::vector<std::vector<float>> Bilerper::GetBidMesh(
    const std::vector<float>& xs, const std::vector<float>& ys) const {
  assert(std::is_sorted(ys.begin(), ys.end()));
  std::vector<std::vector<float>> res(ys.size());
  int start_i = 0;
  int end_i = ys.size() - 1;
  auto boundary = slices_.front()(xs);
  while (start_i <= end_i && ys[start_i] <= y_int_.min) {
    res[start_i] = boundary;
    ++start_i;
  }
  boundary = slices_.back()(xs);
  while (end_i >= start_i && ys[end_i] >= y_int_.max) {
    res[end_i] = boundary;
    --end_i;
  }

  std::vector<float> slice1(xs.size());
  std::vector<float> slice2(xs.size());
  std::vector<float> out_slice(xs.size());
  int old_n = -2;
  for (int i = start_i; i <= end_i; ++i) {
    float n_float = (ys[i] - y_int_.min) / y_spacing_;
    int n = static_cast<int>(n_float);
    float alpha = n_float - n;
    if (n == old_n + 1) {
      slice1 = std::move(slice2);
      slice2 = slices_[n + 1](xs);
    }
    if (n > old_n + 1) {
      slice1 = slices_[n](xs);
      slice2 = slices_[n + 1](xs);
    }
    old_n = n;
    std::transform(
        slice1.begin(), slice1.end(), slice2.begin(), out_slice.begin(),
        [alpha](float y1, float y2) { return (1 - alpha) * y1 + alpha * y2; });
    res[i] = out_slice;
  }
  return res;
}

std::function<float(float)> Bilerper::GenerateSlice(float y) const {
  if (y < y_int_.min) {
    return slices_[0];
  } else if (y >= y_int_.max) {
    return slices_[slices_.size() - 1];
  } else {
    float n_float = (y - y_int_.min) / y_spacing_;
    int n = static_cast<int>(n_float);
    float alpha = n_float - n;
    return [this, alpha, n](float x) {
      return (1 - alpha) * slices_[n](x) + alpha * slices_[n + 1](x);
    };
  }
}

}  // namespace numericaldists
