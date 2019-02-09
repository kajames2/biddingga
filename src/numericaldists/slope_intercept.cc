#include "numericaldists/slope_intercept.h"

#include <vector>

namespace numericaldists {

float EvaluateLine(SlopeIntercept line, float x) {
  return line.slope * x + line.intercept;
}

std::vector<SlopeIntercept> MakeLines(const std::vector<float>& xs,
                                      const std::vector<float>& ys) {
  std::vector<SlopeIntercept> lines;
  lines.reserve(ys.size() - 1);
  for (int i = 1; i < ys.size(); ++i) {
    float slope = (ys[i] - ys[i - 1]) / (xs[i] - xs[i - 1]);
    float intercept = ys[i] - xs[i] * slope;
    lines.push_back({slope, intercept});
  }
  return lines;
}

std::vector<SlopeIntercept> MakeLines(const std::vector<float>& ys,
                                      Interval x_range) {
  std::vector<SlopeIntercept> lines;
  int n_lines = ys.size() - 1;
  lines.reserve(n_lines);
  float x_sep = GetSpan(x_range) / n_lines;
  float x = x_range.min;
  for (int i = 1; i < ys.size(); ++i) {
    x += x_sep;
    float slope = (ys[i] - ys[i - 1]) / x_sep;
    float intercept = ys[i] - x * slope;
    lines.push_back({slope, intercept});
  }
  return lines;
}

}  // namespace numericaldists
