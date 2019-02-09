#ifndef _NUMERICALDISTS_BILERPER_
#define _NUMERICALDISTS_BILERPER_

#include <functional>
#include <vector>

#include "numericaldists/interval.h"
#include "numericaldists/piecewise_linear.h"

namespace numericaldists {

class Bilerper {
 public:
  Bilerper() {}
  Bilerper(Interval x_int, Interval y_int, std::vector<std::vector<float>> zs);
  float operator()(float x, float y) const { return GetBid(x, y); }
  std::vector<std::vector<float>> operator()(
      const std::vector<float>& xs, const std::vector<float>& ys) const {
    return GetBidMesh(xs, ys);
  }
  float GetBid(float x, float y) const;
  std::vector<std::vector<float>> GetBidMesh(
      const std::vector<float>& xs, const std::vector<float>& ys) const;
  std::function<float(float)> GenerateSlice(float y) const;
 private:
  float GetYMax() const;
  int GetIndex(float y) const;
  float GetAlpha(float y) const;
  std::vector<numericaldists::PiecewiseLinear> slices_;
  float y_spacing_;
  Interval y_int_;
};

}  // namespace numericaldists

#endif  // _NUMERICALDISTS_BILERPER_
