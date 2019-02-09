#include "numericaldists/function_ops.h"

#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>
#include <iostream>

#include <boost/math/quadrature/gauss.hpp>

#include "numericaldists/bilerper.h"
#include "numericaldists/line_ops.h"
#include "numericaldists/piecewise_linear.h"
#include "numericaldists/uneven_piecewise_linear.h"

namespace numericaldists {

PiecewiseLinear ResampleFunction(const std::function<float(float)>& func,
                                 Interval interval, int n_samples) {
  auto x_samples = GetMesh(interval, n_samples);
  std::vector<float> y_samples(n_samples);
  std::transform(x_samples.begin(), x_samples.end(), y_samples.begin(), func);
  return PiecewiseLinear(y_samples, interval);
}

PiecewiseLinear ResampleFunction(const UnevenPiecewiseLinear& func,
                                 Interval interval, int n_samples) {
  auto x_samples = GetMesh(interval, n_samples);
  auto y_samples = func(x_samples);
  return PiecewiseLinear(y_samples, interval);
}

UnevenPiecewiseLinear ApproximateInverse(
    const std::function<float(float)>& func, Interval interval, int n_samples) {
  auto x_samples = GetMesh(interval, n_samples);
  std::vector<float> y_samples(n_samples);
  std::transform(x_samples.begin(), x_samples.end(), y_samples.begin(), func);
  if (y_samples.front() > y_samples.back()) {
    std::reverse(x_samples.begin(), x_samples.end());
    std::reverse(y_samples.begin(), y_samples.end());
  }
  return UnevenPiecewiseLinear(y_samples, x_samples);
}

PiecewiseLinear ApproximateDerivative(const std::function<float(float)>& func,
                                      Interval interval, int n_samples) {
  float d = GetSpan(interval) / n_samples;
  auto xd_samples = GetMesh(interval, n_samples + 1);
  auto x_samples =
      GetMesh(Interval{interval.min + d / 2, interval.max - d / 2}, n_samples);
  std::vector<float> yd_samples(n_samples + 1);
  std::transform(xd_samples.begin(), xd_samples.end(), yd_samples.begin(),
                 func);
  std::vector<float> y_samples(n_samples);
  for (int i = 0; i < n_samples; ++i) {
    y_samples[i] = (yd_samples[i + 1] - yd_samples[i]) / d;
  }
  return PiecewiseLinear(y_samples, {x_samples.front(), x_samples.back()});
}

// Returns areas under curve for between each samples.  out[0] = 0, and out[i] =
// area between x_sample[i-1] and x_sample[i]
std::vector<float> ApproximateAreas(
    std::function<float(float)> func, const std::vector<float>& x_samples,
    std::function<float(float)> integrand_func) {
  std::vector<float> int_samples(x_samples.size());
  int_samples[0] = 0;
  for (int i = 1; i < int_samples.size(); ++i) {
    int_samples[i] = boost::math::quadrature::gauss<float, 7>::integrate(
        [func, integrand_func](float x) { return integrand_func(x) * func(x); },
        x_samples[i - 1], x_samples[i]);
  }
  return int_samples;
}

// Returns a functor, f, for which over the interval [a,b), f(x) will return
// int_a^x{integrand_func * func(x)dx}
PiecewiseLinear ApproximateIntegralBelow(
    std::function<float(float)> func, Interval interval, int n_samples,
    std::function<float(float)> integrand_func) {
  auto x_samples = GetMesh(interval, n_samples);
  auto area_samples = ApproximateAreas(func, x_samples, integrand_func);
  std::partial_sum(area_samples.begin(), area_samples.end(),
                   area_samples.begin());
  return PiecewiseLinear(area_samples, {x_samples.front(), x_samples.back()});
}

// Returns a functor, f, for which over the interval [a,b), f(x) will return
// int_x^b{integrand_func * func(x)dx}
PiecewiseLinear ApproximateIntegralAbove(
    std::function<float(float)> func, Interval interval, int n_samples,
    std::function<float(float)> integrand_func) {
  auto x_samples = GetMesh(interval, n_samples);
  auto area_samples = ApproximateAreas(func, x_samples, integrand_func);
  std::rotate(area_samples.begin(), area_samples.begin() + 1,
              area_samples.end());
  std::partial_sum(area_samples.rbegin(), area_samples.rend(),
                   area_samples.rbegin());
  return PiecewiseLinear(area_samples, {x_samples.front(), x_samples.back()});
}

Bilerper ResampleFunction2D(const std::function<float(float, float)>& func,
                            Interval x_int, Interval y_int, int n_samples) {
  auto x_samples = GetMesh(x_int, n_samples);
  auto y_samples = GetMesh(y_int, n_samples);
  std::vector<std::vector<float>> z_slices;
  for (auto y_it = y_samples.rbegin(); y_it != y_samples.rend(); ++y_it) {
    float y = *y_it;
    std::vector<float> z_samples(n_samples);
    std::transform(x_samples.begin(), x_samples.end(), z_samples.begin(),
                   [y, &func](float x) { return func(x, y); });
    z_slices.push_back(z_samples);
  }
  return Bilerper(x_int, y_int, z_slices);
}

}  // namespace numericaldists
