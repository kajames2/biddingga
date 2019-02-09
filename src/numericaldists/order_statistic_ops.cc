#include "numericaldists/order_statistic_ops.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

#include "numericaldists/bilerper.h"
#include "numericaldists/combination_generation.h"
#include "numericaldists/function_ops.h"
#include "numericaldists/interval.h"
#include "numericaldists/line_ops.h"
#include "numericaldists/piecewise_linear.h"

#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>

namespace numericaldists {

std::vector<float> ApproximateKthLowestOrderStatisticCDFHelper(
    const std::vector<std::function<float(float)>>& cdfs, Interval interval,
    int k, int n_samples) {
  int n = cdfs.size();
  assert(n >= 1 && k <= n);
  auto x_samples = GetMesh(interval, n_samples);
  std::vector<float> cdf_samples(n_samples, 0);

  auto comb = GetFirstCanonicalCombination(k);
  while (comb < (1 << n)) {
    for (int i = 0; i < n_samples; ++i) {
      float prob = 1;
      for (int j = 0; j < n; ++j) {
        float prob_below = cdfs[j](x_samples[i]);
        unsigned int is_below = (comb >> j) & 1;
        prob *= is_below * prob_below + (1 - is_below) * (1 - prob_below);
      }
      cdf_samples[i] += prob;
    }
    comb = GetNextCanonicalCombination(comb);
  }

  if (k < n) {
    auto rest_cdf = ApproximateKthLowestOrderStatisticCDFHelper(
        cdfs, interval, k + 1, n_samples);
    std::transform(cdf_samples.begin(), cdf_samples.end(), rest_cdf.begin(),
                   cdf_samples.begin(), std::plus<float>());
  }
  return cdf_samples;
}

PiecewiseLinear ApproximateKthLowestOrderStatisticCDF(
    const std::vector<std::function<float(float)>>& cdfs, Interval interval,
    int k, int n_samples) {
  auto cdf_samples =
      ApproximateKthLowestOrderStatisticCDFHelper(cdfs, interval, k, n_samples);
  return PiecewiseLinear(cdf_samples, {interval.min, interval.max});
}

PiecewiseLinear ApproximateKthLowestOrderStatisticPDF(
    const std::vector<std::function<float(float)>>& cdfs, Interval interval,
    int k, int n_samples) {
  auto cdf =
      ApproximateKthLowestOrderStatisticCDF(cdfs, interval, k, n_samples);
  return ApproximateDerivative(cdf, interval);
}

PiecewiseLinear ApproximateHighestOrderStatisticCDF(
    const std::vector<std::function<float(float)>>& cdfs, Interval interval,
    int n_samples) {
  return ApproximateKthLowestOrderStatisticCDF(cdfs, interval, cdfs.size(),
                                               n_samples);
}

PiecewiseLinear ApproximateKthLowestOrderStatisticCDF(
    const std::function<float(float)>& cdf, Interval interval, int n_draws,
    int k, int n_samples) {
  auto pdf = ApproximateKthLowestOrderStatisticPDF(cdf, interval, n_draws, k,
                                                   n_samples);
  return ApproximateIntegralBelow(pdf, interval, n_samples);
}

PiecewiseLinear ApproximateKthLowestOrderStatisticPDF(
    const std::function<float(float)>& cdf, Interval interval, int n_draws,
    int k, int n_samples) {
  assert(n_draws > 0 && k > 0 && k <= n_draws);
  auto pdf = ApproximateDerivative(cdf, interval);
  auto x_samples = GetMesh(interval, n_samples);
  float bin = boost::math::binomial_coefficient<float>(n_draws, k);
  std::vector<float> pdf_samples(n_samples);
  std::transform(x_samples.begin(), x_samples.end(), pdf_samples.begin(),
                 [bin, k, &cdf, &pdf, n_draws](float x) {
                   return bin * k * std::pow(cdf(x), k - 1) *
                          std::pow(1 - cdf(x), n_draws - k) * pdf(x);
                 });
  return PiecewiseLinear(pdf_samples, {x_samples.front(), x_samples.back()});
}

PiecewiseLinear ApproximateHighestOrderStatisticCDF(
    const std::function<float(float)> cdf, Interval interval, int n_draws,
    int n_samples) {
  return ApproximateKthLowestOrderStatisticCDF(cdf, interval, n_draws, n_draws,
                                               n_samples);
}

Bilerper ApproximateJointOrderStatisticPDF(const Distribution& dist,
                                           int n_draws, int j, int k,
                                           int n_samples) {
  assert(n_draws > 1);
  assert(k > j);
  Interval interval{lower(dist), upper(dist)};
  auto x_mesh = GetMesh(interval, n_samples);
  auto y_mesh = x_mesh;
  std::vector<std::vector<float>> zs;
  auto fact = boost::math::factorial<double>;
  double n_combs =
      fact(n_draws) / (fact(j - 1) * fact(k - j - 1) * fact(n_draws - k));
  for (auto y_ptr = y_mesh.rbegin(); y_ptr != y_mesh.rend(); ++y_ptr) {
    auto y = *y_ptr;
    std::vector<float> z_slice(n_samples);
    std::transform(
        x_mesh.begin(), x_mesh.end(), z_slice.begin(),
        [&dist, n_draws, y, j, k, n_combs](float x) {
          if (x > y) {
            return 0.;
          }
          float Fx = cdf(dist, x);
          float Fy = cdf(dist, y);
          return n_combs * std::pow(Fx, j - 1) * std::pow(Fy - Fx, k - 1 - j) *
                 std::pow(1 - Fx, n_draws - k) * pdf(dist, x) * pdf(dist, y);
        });
    zs.push_back(z_slice);
  }
  return Bilerper(interval, interval, zs);
}

Bilerper ApproximateLowestHighestJointOrderStatisticPDF(
    const Distribution& dist, int n_draws, int n_samples) {
  return ApproximateJointOrderStatisticPDF(dist, n_draws, 1, n_draws,
                                           n_samples);
}

Bilerper ApproximateLowestHighestJointOrderStatisticCDF(
    const std::function<float(float)> cdf, Interval interval, int n_draws,
    int n_samples) {
  assert(n_draws > 1);
  auto x_mesh = GetMesh(interval, n_samples);
  auto y_mesh = x_mesh;
  std::vector<std::vector<float>> zs;
  for (auto y_ptr = y_mesh.rbegin(); y_ptr != y_mesh.rend(); ++y_ptr) {
    float y = *y_ptr;
    std::vector<float> z_slice(n_samples);
    std::transform(x_mesh.begin(), x_mesh.end(), z_slice.begin(),
                   [&cdf, n_draws, y](float x) {
                     return std::pow(cdf(y), n_draws) -
                            std::pow(cdf(y) - cdf(x), n_draws);
                   });
    zs.push_back(z_slice);
  }
  return Bilerper(interval, interval, zs);
}

}  // namespace numericaldists
