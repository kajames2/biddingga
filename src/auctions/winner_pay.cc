#include "auctions/winner_pay.h"

#include <algorithm>
#include <iostream>

#include <boost/math/distributions/uniform.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include "numericaldists/distribution_ops.h"
#include "numericaldists/function_ops.h"
#include "numericaldists/order_statistic_ops.h"

using namespace boost::math;
using namespace numericaldists;

namespace auctions {

WinnerPay::WinnerPay(std::vector<Distribution> value_dists, int pay_nth_highest)
    : bid_funcs_(value_dists.size()),
      others_highest_order_stat_(value_dists.size()),
      others_nth_highest_exp_bid_(value_dists.size()),
      value_dists_(value_dists),
      n_players_(value_dists.size()),
      pre_calculated_(false),
      pay_nth_highest_(pay_nth_highest),
      max_bid_(upper(value_dists)) {}

WinnerPay::WinnerPay(std::vector<Distribution> value_dists, int pay_nth_highest,
                     float max_bid)
    : bid_funcs_(value_dists.size()),
      others_highest_order_stat_(value_dists.size()),
      others_nth_highest_exp_bid_(value_dists.size()),
      value_dists_(value_dists),
      n_players_(value_dists.size()),
      pre_calculated_(false),
      pay_nth_highest_(pay_nth_highest),
      max_bid_(max_bid) {}

void WinnerPay::AcceptStrategy(PiecewiseLinear bid_func, int id) {
  bid_funcs_[id] = std::move(bid_func);
  pre_calculated_ = false;
}

float WinnerPay::GetFitness(const PiecewiseLinear& bid_func, int id) const {
  if (!pre_calculated_) {
    Precalculate();
  }

  float exp_profit = quadrature::gauss_kronrod<float, 121>::integrate(
      [this, &bid_func, id](float value) {
        return GetIntegrand(bid_func, id, value);
      },
      lower(value_dists_[id]), upper(value_dists_[id]), 0, 0);
  return exp_profit;
}

void WinnerPay::Precalculate() const {
  Interval interval = Interval{0, max_bid_};
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
    others_highest_order_stat_[i] =
        ApproximateHighestOrderStatisticCDF(other_cdfs, interval);
    if (pay_nth_highest_ > 1) {
      auto nth_order_stat = ApproximateKthLowestOrderStatisticCDF(
          other_cdfs, interval, n_players_ - pay_nth_highest_ + 1);
      others_nth_highest_exp_bid_[i] =
          ApproximateCDFExpectedValueFunction(nth_order_stat, interval);
    }
  }
  pre_calculated_ = true;
}

double WinnerPay::GetIntegrand(const PiecewiseLinear& bid_func, int id,
                               float value) const {
  double bid = bid_func(value);
  double prob_value = pdf(value_dists_[id], value);
  double prob_win_or_tie = others_highest_order_stat_[id](bid);
  double prob_win_outright =
      others_highest_order_stat_[id](bid - max_bid_ / 1000);
  double prob_win = (prob_win_or_tie + prob_win_outright) / 2;
  if (pay_nth_highest_ == 1) {
    return (value - bid) * prob_win * prob_value;
  } else {
    double exp_payment_given_win = others_nth_highest_exp_bid_[id](bid);
    return (value - exp_payment_given_win) * prob_win * prob_value;
  }
}

}  // namespace auctions
