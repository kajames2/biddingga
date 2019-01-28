#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>

#include "numericaldists/distribution.h"
#include "numericaldists/function_ops.h"
#include "numericaldists/line_ops.h"
#include "numericaldists/piecewise_linear.h"
#include "numericaldists/uneven_piecewise_linear.h"

namespace numericaldists {

UnevenPiecewiseLinear ApproximateRandomVariableFunctionCDF(
    const Distribution& dist, const std::function<float(float)>& func,
    int n_samples) {
  auto prob_samples = GetMesh(Interval{0.000001, 0.999999}, n_samples);
  std::vector<float> quant_samples(n_samples);
  std::transform(prob_samples.begin(), prob_samples.end(),
                 quant_samples.begin(),
                 [&dist](float x) { return quantile(dist, x); });
  std::vector<float> f_samples(n_samples);
  std::transform(quant_samples.begin(), quant_samples.end(), f_samples.begin(),
                 func);
  std::sort(f_samples.begin(), f_samples.end());
  auto segments = PointsToLines(f_samples, prob_samples);
  return UnevenPiecewiseLinear(segments);
}

PiecewiseLinear ApproximateExpectedValueFunction(
    std::function<float(float)> pdf, std::function<float(float)> cdf,
    Interval interval, int n_samples) {
  auto x_samples = GetMesh(interval, n_samples);
  auto area_samples =
      ApproximateAreas(pdf, x_samples, [](float x) { return x; });
  std::partial_sum(area_samples.begin(), area_samples.end(),
                   area_samples.begin());
  area_samples[0] = x_samples[0];
  for (int i = 1; i < n_samples; ++i) {
    area_samples[i] =
        std::min(area_samples[i] / cdf(x_samples[i]), x_samples[i]);
  }
  auto segments = PointsToLines(x_samples, area_samples);
  return PiecewiseLinear(segments);
}

PiecewiseLinear ApproximatePDFExpectedValueFunction(
    std::function<float(float)> pdf, Interval interval, int n_samples) {
  auto cdf = ApproximateIntegralBelow(pdf, interval, n_samples);
  return ApproximateExpectedValueFunction(pdf, cdf, interval, n_samples);
}

PiecewiseLinear ApproximateCDFExpectedValueFunction(
    std::function<float(float)> cdf, Interval interval, int n_samples) {
  auto pdf = ApproximateDerivative(cdf, interval, n_samples);
  return ApproximateExpectedValueFunction(pdf, cdf, interval, n_samples);
}

PiecewiseLinear ApproximateConditionalXPDF(
    const std::function<float(float, float)>& pdf, Interval x_int, float y,
    int n_samples) {
  auto x_mesh = GetMesh(x_int, n_samples);
  std::vector<float> cond(n_samples, 0);
  std::transform(x_mesh.begin(), x_mesh.end(), cond.begin(),
                 [y, &pdf](float x) { return pdf(x, y); });
  auto total = std::accumulate(cond.begin(), cond.end(), 0);
  for (auto& c : cond) c /= total;
  return PiecewiseLinear(x_mesh, cond);
}

PiecewiseLinear ApproximateConditionalYPDF(
    const std::function<float(float, float)>& pdf, Interval y_int, float x,
    int n_samples) {
  auto y_mesh = GetMesh(y_int, n_samples);
  std::vector<float> cond(n_samples, 0);
  std::transform(y_mesh.begin(), y_mesh.end(), cond.begin(),
                 [x, &pdf](float y) { return pdf(x, y); });
  auto total = std::accumulate(cond.begin(), cond.end(), 0);
  for (auto& c : cond) c /= total;
  return PiecewiseLinear(y_mesh, cond);
}

PiecewiseLinear ApproximateRandomVariableFunctionCDF(
    const std::function<float(float, float)>& pdf,
    const std::function<float(float, float)>& func, Interval x_int,
    Interval y_int, Interval f_int, int n_samples) {
  float x_size = GetSpan(x_int) / n_samples;
  float y_size = GetSpan(y_int) / n_samples;
  float section_area = x_size * y_size;
  auto x_mesh = GetMesh(
      Interval{x_int.min + x_size / 2, x_int.max - x_size / 2}, n_samples);
  auto y_mesh = GetMesh(
      Interval{y_int.min + y_size / 2, y_int.max - y_size / 2}, n_samples);
  std::vector<float> f_mesh(n_samples, 0);
  float f_bin_size = GetSpan(f_int) / n_samples;
  for (auto x : x_mesh) {
    for (auto y : y_mesh) {
      float f_val = func(x, y);
      int bin = static_cast<int>((f_val - f_int.min) / f_bin_size);
      assert(bin >= 0 && bin < n_samples);
      f_mesh[bin] += pdf(x, y) * section_area;
    }
  }
  std::partial_sum(f_mesh.begin(), f_mesh.end(), f_mesh.begin());
  return PiecewiseLinear(f_mesh, f_int);
}

// Must be monotonic transformations
std::function<float(float, float)> MultivariatePDFDomainTransform(
    const std::function<float(float, float)>& func,
    const std::function<float(float, float)>& x_trans,
    const std::function<float(float, float)>& y_trans,
    const std::function<float(float, float)>& jacobian, Interval x_int,
    Interval y_int) {
  return [func, x_trans, y_trans, jacobian, x_int, y_int](float u, float v) -> float{
    float x = x_trans(u,v);
    float y = y_trans(u,v);
    if (!InInterval(x_int, x) || !InInterval(y_int, y)) {return 0;}
    return func(x_trans(u, v), y_trans(u, v)) * jacobian(u, v);
  };
}

}  // namespace numericaldists
