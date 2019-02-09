#ifndef _NUMERICALDISTS_SLOPE_INTERCEPT_H_
#define _NUMERICALDISTS_SLOPE_INTERCEPT_H_

#include <vector>

#include "numericaldists/interval.h"

namespace numericaldists {

struct SlopeIntercept {
  float slope;
  float intercept;
};

float EvaluateLine(SlopeIntercept line, float x);
std::vector<SlopeIntercept> MakeLines(const std::vector<float>& xs,
                                      const std::vector<float>& ys);
std::vector<SlopeIntercept> MakeLines(const std::vector<float>& ys,
                                      Interval x_range);

}  // namespace numericaldists

#endif  // _NUMERICALDISTS_SLOPE_INTERCEPT_H_
