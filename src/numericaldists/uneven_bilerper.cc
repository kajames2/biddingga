#include "numericaldists/uneven_bilerper.h"

#include <cmath>
#include <functional>
#include <vector>

#include "numericaldists/interval.h"
#include "numericaldists/uneven_piecewise_linear.h"

namespace numericaldists {

UnevenBilerper::UnevenBilerper(std::vector<float> xs, std::vector<float> ys,
                               std::vector<std::vector<float>> zs)
    : ys_(ys) {
  for (auto slice = zs.rbegin(); slice != zs.rend(); ++slice) {
    slices_.push_back(UnevenPiecewiseLinear(xs, *slice));
  }
}

UnevenBilerper::UnevenBilerper(std::vector<float> ys,
                               std::vector<std::function<float(float)>> slices)
    : slices_(slices), ys_(ys) {}

float UnevenBilerper::operator()(float x, float y) const {
  if (y <= ys_[0]) {
    return slices_[0](x);
  } else if (y >= ys_.back()) {
    return slices_[slices_.size() - 1](x);
  } else {
    int index = GetIndex(y);
    float alpha = GetAlpha(index, y);
    return (1 - alpha) * slices_[index - 1](x) + alpha * slices_[index](x);
  }
}

int UnevenBilerper::GetIndex(float y) const {
  for (int i = 1; i < ys_.size(); ++i) {
    if (y <= ys_[i]) {
      return i;
    }
  }
  return ys_.size() - 1;
}

float UnevenBilerper::GetAlpha(int i, float y) const {
  return (y - ys_[i - 1]) / (ys_[i] - ys_[i - 1]);
}

}  // namespace numericaldists
