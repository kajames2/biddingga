#include "auctions/second_price.h"

#include <algorithm>

#include <boost/math/distributions/uniform.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include "numericaldists/function_ops.h"
#include "numericaldists/order_statistic_ops.h"
#include "numericaldists/distribution_ops.h"

using namespace boost::math;
using namespace numericaldists;

namespace auctions {

SecondPrice::SecondPrice(std::vector<Distribution> value_dists)
    : bid_funcs_(value_dists.size()),
      order_stat_funcs_(value_dists.size()),
      exp_value_funcs_(value_dists.size()),
      value_dists_(value_dists),
      n_players_(value_dists.size()),
      pre_calculated_(false) {}

void SecondPrice::AcceptStrategy(std::function<float(float)> bid_func,
                                        int id) {
  bid_funcs_[id] = bid_func;
  pre_calculated_ = false;
}

float SecondPrice::GetFitness(
    const std::function<float(float)>& bid_func, int id) const {
  if (!pre_calculated_) {
    Precalculate();
  }
  float exp_profit = quadrature::gauss_kronrod<float, 61>::integrate(
      [this, &bid_func, id](float value) {
        return GetIntegrand(bid_func, id, value);
      },
      lower(value_dists_[id]), upper(value_dists_[id]), 0, 0);
  return exp_profit;
}

void SecondPrice::Precalculate() const {
  Interval interval = Interval{0, upper(value_dists_)};
  std::vector<std::function<float(float)>> cdfs;
  for (int i = 0; i < n_players_; ++i) {
    cdfs.push_back(ResampleFunction(
        ApproximateRandomVariableFunctionCDF(value_dists_[i], bid_funcs_[i]),
        interval));
  }

  for (int i = 0; i < n_players_; ++i) {
    std::vector<std::function<float(float)>> other_cdfs;
    std::vector<Distribution> other_dists;
    for (int j = 0; j < n_players_; ++j) {
      if (i != j) {
        other_cdfs.emplace_back(cdfs[j]);
        other_dists.emplace_back(value_dists_[j]);
      }
    }
    order_stat_funcs_[i] = ApproximateKthLowestOrderStatisticCDF(
        other_cdfs, interval, n_players_ - 1);
    exp_value_funcs_[i] =
        ApproximateCDFExpectedValueFunction(order_stat_funcs_[i], interval);
  }
  pre_calculated_ = true;
}

double SecondPrice::GetIntegrand(
    const std::function<float(float)>& bid_func, int id, float value) const {
  double bid = bid_func(value);
  double prob_bid = pdf(value_dists_[id], value);
  double prob_win = order_stat_funcs_[id](bid);
  double exp_second_bid_given_win = exp_value_funcs_[id](bid);
  return (value - exp_second_bid_given_win) * prob_win * prob_bid;
}

}  // namespace auctions
