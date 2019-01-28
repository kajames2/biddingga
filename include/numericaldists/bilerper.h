#ifndef _NUMERICALDISTS_BILERPER_
#define _NUMERICALDISTS_BILERPER_

#include <functional>
#include <vector>

#include "numericaldists/interval.h"

namespace numericaldists {

class Bilerper {
 public:
  Bilerper() {}
  Bilerper(Interval x_int, Interval y_int, std::vector<std::vector<float>> zs);
  Bilerper(Interval x_int,
           std::vector<std::function<float(float)>> slices);
  float operator()(float x, float y) const;

 private:
  float GetYMax() const;
  int GetIndex(float y) const;
  float GetAlpha(float y) const;
  std::vector<std::function<float(float)>> slices_;
  float y_spacing_;
  float y_min_;
};

}  // namespace numericaldists

#endif  // _NUMERICALDISTS_BILERPER_
