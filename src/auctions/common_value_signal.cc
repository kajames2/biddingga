#include "auctions/common_value_signal.h"

#include <iostream>

#include <algorithm>

#include <boost/math/distributions/uniform.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include "numericaldists/distribution_ops.h"
#include "numericaldists/function_ops.h"
#include "numericaldists/order_statistic_ops.h"

using namespace boost::math;
using namespace numericaldists;

namespace auctions {

CommonValueSignal::CommonValueSignal(std::vector<int> n_draws, float epsilon,
                                     Interval discount_int)
    : n_players_(n_draws.size()),
      pre_calculated_(false),
      n_draws_(n_draws),
      epsilon_(epsilon),
      discount_funcs_(n_draws.size()),
      others_bids_cdfs_(n_draws.size()),
      signal_int_(Interval{-epsilon, epsilon}),
      range_int_(Interval{0, 2 * epsilon}),
      bid_int_(
          Interval{-epsilon - discount_int.max, epsilon - discount_int.min}) {
  value_dists_ =
      std::vector<std::function<float(float, float)>>(n_draws.size());
  for (int i = 0; i < n_players_; ++i) {
    int draw = n_draws_[i];
    if (draw > 1) {
      auto joint = ApproximateJointOrderStatisticPDF(
          uniform_distribution<>(signal_int_.min, signal_int_.max), draw, 1,
          draw);
      auto mstar_range = MultivariatePDFDomainTransform(
          joint,
          [epsilon](float mstar, float range) {
            return -epsilon + mstar + 0.5 * range;
          },
          [epsilon](float mstar, float range) {
            return epsilon + mstar - 0.5 * range;
          },
          [](float mstar, float range) { return 1; }, signal_int_, signal_int_);
      auto mstar_range_resampled =
          ResampleFunction2D(mstar_range, signal_int_, range_int_);
      value_dists_[i] = mstar_range_resampled;
    } else {
      float inverse_max_range = 1 / (2 * epsilon_);
      value_dists_[i] = [inverse_max_range](float mstar, float range) {
        return inverse_max_range;
      };
    }
  }
}

void CommonValueSignal::AcceptStrategy(
    std::function<float(float)> discount_func, int id) {
  discount_funcs_[id] = discount_func;
  pre_calculated_ = false;
}

float CommonValueSignal::GetFitness(
    const std::function<float(float)>& discount_func, int id) const {
  if (!pre_calculated_) {
    Precalculate();
  }

  float exp_profit = 0;
  if (n_draws_[id] > 1) {
    exp_profit = quadrature::gauss_kronrod<float, 61>::integrate(
        [this, &discount_func, id](float range) {
          return GetIntegrand(discount_func, id, range);
        },
        range_int_.min, range_int_.max, 0, 0);
  } else {
    exp_profit = GetIntegrand(discount_func, id, 2 * epsilon_);
  }
  return exp_profit;
}

double CommonValueSignal::GetIntegrand(
    const std::function<float(float)>& discount_func, int id,
    float range) const {
  return quadrature::gauss_kronrod<float, 61>::integrate(
      [this, range, &discount_func, id](float m_star) {
        float bid = m_star - discount_func(range);
        return (0 - bid) * others_bids_cdfs_[id](bid) *
               value_dists_[id](m_star, range);
      },
      signal_int_.min, signal_int_.max, 0, 0);
}

void CommonValueSignal::Precalculate() const {
  std::vector<std::function<float(float)>> cdfs;
  for (int i = 0; i < n_players_; ++i) {
    if (n_draws_[i] > 1) {
      std::function<float(float)> cdf = ApproximateRandomVariableFunctionCDF(
          value_dists_[i],
          [this, i](float m_star, float range) {
            return m_star - discount_funcs_[i](range);
          },
          signal_int_, range_int_, bid_int_);
      cdfs.push_back(ResampleFunction(cdf, bid_int_));
    } else {
      std::function<float(float)> cdf = ApproximateRandomVariableFunctionCDF(
          uniform_distribution<>(signal_int_.min, signal_int_.max),
          [this, i](float range) { return 0 - discount_funcs_[i](range); });
    }
  }

  for (int i = 0; i < n_players_; ++i) {
    std::vector<std::function<float(float)>> other_cdfs;
    for (int j = 0; j < n_players_; ++j) {
      if (i != j) {
        other_cdfs.emplace_back(cdfs[j]);
      }
    }
    others_bids_cdfs_[i] = ApproximateKthLowestOrderStatisticCDF(
        other_cdfs, bid_int_, n_players_ - 1);
  }
  pre_calculated_ = true;
}

}  // namespace auctions
