#include "numericaldists/bilerper.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

#include "numericaldists/piecewise_linear.h"

namespace numericaldists {

Bilerper::Bilerper(Interval x_int, Interval y_int,
                   std::vector<std::vector<float>> zs)
    : y_spacing_(GetSpan(y_int)/(zs.size() - 1)), y_min_(y_int.min) {
  for (auto slice = zs.rbegin(); slice != zs.rend(); ++slice) {
    slices_.push_back(PiecewiseLinear(*slice, x_int));
  }
}

Bilerper::Bilerper(Interval y_int,
                   std::vector<std::function<float(float)>> slices)
    : slices_(slices), y_spacing_(GetSpan(y_int)/(slices.size() - 1)), y_min_(y_int.min) {
  std::reverse(slices_.begin(), slices_.end());
}

float Bilerper::operator()(float x, float y) const {
  if (y < y_min_) {
    return slices_[0](x);
  } else if (y >= GetYMax()) {
    return slices_[slices_.size() - 1](x);
  } else {
    int index = GetIndex(y);
    float alpha = (y - y_min_ - (y_spacing_ * index)) / y_spacing_;
    //float alpha = GetAlpha(y);
    return (1 - alpha) * slices_[index](x) + alpha * slices_[index + 1](x);
  }
}

float Bilerper::GetYMax() const {
  return y_min_ + (slices_.size() - 1) * y_spacing_;
}
int Bilerper::GetIndex(float y) const {
  return static_cast<int>((y - y_min_) / y_spacing_);
}
float Bilerper::GetAlpha(float y) const {
  return std::fmod(y - y_min_, y_spacing_) / y_spacing_;
}

}  // namespace numericaldists
