#include "numericaldists/order_statistic_ops.h"

#include <cassert>
#include <cmath>
#include <vector>
#include <iostream>

#include "numericaldists/combination_generation.h"
#include "numericaldists/function_ops.h"

#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <eigen3/Eigen/Core>

using namespace Eigen;
namespace numericaldists {

ArrayXd KthLowestOrderStatisticCDF(const ArrayXd& xs,
                                   const std::vector<ArrayXd>& cdfs, int k) {
  int n = cdfs.size();
  assert(n >= 1 && k <= n);
  int size = cdfs[0].size();
  ArrayXd cdf_samples = ArrayXd::Zero(size);
  auto comb = GetFirstCanonicalCombination(k);
  while (comb < (1 << n)) {
    ArrayXd probs = ArrayXd::Ones(size);
    for (int j = 0; j < n; ++j) {
      unsigned int is_below = (comb >> j) & 1;
      if (is_below == 1) {
        probs *= cdfs[j];
      } else {
        probs *= (1 - cdfs[j]);
      }
    }
    cdf_samples += probs;
    comb = GetNextCanonicalCombination(comb);
  }
  if (k < n) {
    cdf_samples += KthLowestOrderStatisticCDF(xs, cdfs, k + 1);
  }
  return cdf_samples;
}

ArrayXd KthLowestOrderStatisticPDF(const ArrayXd& xs,
                                   const std::vector<ArrayXd>& cdfs, int k) {
  ArrayXd cdf = KthLowestOrderStatisticCDF(xs, cdfs, k);
  return Derivative(xs, cdf);
}

ArrayXd HighestOrderStatisticCDF(const ArrayXd& xs,
                                 const std::vector<ArrayXd>& cdfs) {
  return KthLowestOrderStatisticCDF(xs, cdfs, cdfs.size());
}

ArrayXd KthLowestOrderStatisticPDF(const ArrayXd& xs, const ArrayXd& cdf,
                                   int n_draws, int k) {
  assert(n_draws > 0 && k > 0 && k <= n_draws);
  auto pdf = Derivative(xs, cdf);
  double bin = boost::math::binomial_coefficient<double>(n_draws, k);
  ArrayXd pdf_samples =
      bin * k * cdf.pow(k - 1) * (1 - cdf).pow(n_draws - k) * pdf;
  return pdf_samples;
}

ArrayXd KthLowestOrderStatisticCDF(const ArrayXd& xs, const ArrayXd& cdf,
                                   int n_draws, int k) {
  auto pdf = KthLowestOrderStatisticPDF(xs, cdf, n_draws, k);
  return IntegralBelow(xs, pdf);
}

ArrayXd HighestOrderStatisticCDF(const ArrayXd& xs, const ArrayXd& cdf,
                                 int n_draws) {
  return KthLowestOrderStatisticCDF(xs, cdf, n_draws, n_draws);
}

ArrayXXd JointOrderStatisticPDF(const ArrayXd& xs, const ArrayXd& cdf,
                                int n_draws, int j, int k) {
  assert(n_draws > 1);
  assert(k > j);
  ArrayXd pdf = Derivative(xs, cdf);
  auto fact = boost::math::factorial<double>;
  double n_combs =
      fact(n_draws) / (fact(j - 1) * fact(k - j - 1) * fact(n_draws - k));
  ArrayXd y_ind_pdf =
      n_combs * cdf.pow(j - 1) * (1 - cdf).pow(n_draws - k) * pdf;

  ArrayXXd zs(xs.size(), xs.size());
  for (int i = 0; i < xs.size(); ++i) {
    double y = xs(i);
    double Fy = cdf(i);
    double fy = pdf(i);
    zs.row(i) =
        (xs < y).cast<double>() * y_ind_pdf * (Fy - cdf).pow(k - 1 - j) * fy;
  }
  return zs;
}

ArrayXXd LowestHighestJointOrderStatisticPDF(const ArrayXd& xs,
                                             const ArrayXd& cdf, int n_draws) {
  return JointOrderStatisticPDF(xs, cdf, n_draws, 1, n_draws);
}

ArrayXXd LowestHighestJointOrderStatisticCDF(const ArrayXd& xs,
                                             const ArrayXd& cdf, int n_draws) {
  assert(n_draws > 1);
  ArrayXXd zs(xs.size(), xs.size());
  for (int i = 0; i < xs.size(); ++i) {
    zs.row(i) = std::pow(cdf(i), n_draws) - (cdf(i) - cdf).pow(n_draws);
  }
  return zs;
}

// Joint order statistic for N different distributions...

}  // namespace numericaldists
