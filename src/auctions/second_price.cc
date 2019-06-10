#include "auctions/second_price.h"

#include <algorithm>
#include <iostream>

#include "numericaldists/distribution_ops.h"
#include "numericaldists/function_ops.h"
#include "numericaldists/order_statistic_ops.h"

using namespace numericaldists;
using namespace Eigen;

namespace auctions {

SecondPrice::SecondPrice(std::vector<Distribution> value_dists,
                         int n_internal_samples)
    : bid_cdfs_(value_dists.size()),
      other_highest_cdfs_(value_dists.size()),
      exp_value_funcs_(value_dists.size()),
      n_players_(value_dists.size()),
      n_internal_samples_(n_internal_samples),
      pre_calculated_(false) {
  utility_funcs_ = std::vector<std::function<double(double)>>(
      n_players_, [](double val) { return val; });
  prob_weight_funcs_ = std::vector<std::function<double(double)>>(
      n_players_, [](double prob) { return prob; });
  internal_values_ = ArrayXd::LinSpaced(n_internal_samples, lower(value_dists),
                                        upper(value_dists));
  for (const auto& dist : value_dists) {
    ArrayXd pdfs = internal_values_.unaryExpr(
        [&dist](double x) -> double { return pdf(dist, x); });
    value_pdfs_.push_back(pdfs);
  }
  internal_bids_ =
      ArrayXd::LinSpaced(n_internal_samples_, 0, upper(value_dists));
}

void SecondPrice::AcceptStrategy(const Scatter& bids, int id) {
  ArrayXd interp_bids = Interpolate(bids, internal_values_);
  bid_cdfs_[id] = RandomVariableFunctionCDF(internal_values_, value_pdfs_[id],
                                            interp_bids, internal_bids_);
  pre_calculated_ = false;
}

float SecondPrice::GetFitness(const Scatter& bids_in, int id) const {
  if (!pre_calculated_) {
    Precalculate();
  }

  ArrayXd integrate_vals =
      ArrayXd::LinSpaced((bids_in.xs.size() - 1) * 3, bids_in.xs(0),
                         bids_in.xs(bids_in.xs.size() - 1));
  ArrayXd bids = Interpolate(bids_in, integrate_vals);
  ArrayXd win_probs =
      Interpolate(internal_bids_, other_highest_cdfs_[id], bids);
  ArrayXd exp_second_bid_given_win =
      Interpolate(internal_bids_, exp_value_funcs_[id], bids);
  ArrayXd profits = integrate_vals - exp_second_bid_given_win;
  ArrayXd utils =
      profits.unaryExpr(utility_funcs_[id]) *
          win_probs.unaryExpr(prob_weight_funcs_[id]) +
      utility_funcs_[id](0) * (1 - win_probs).unaryExpr(prob_weight_funcs_[id]);
  ArrayXd likelihoods =
      Interpolate(internal_values_, value_pdfs_[id], integrate_vals);
  return Areas(integrate_vals, utils * likelihoods).sum();
}

float SecondPrice::GetRevenue(const Scatter& bids_in, int id) const {
  if (!pre_calculated_) {
    Precalculate();
  }

  ArrayXd integrate_vals =
      ArrayXd::LinSpaced((bids_in.xs.size() - 1) * 3, bids_in.xs(0),
                         bids_in.xs(bids_in.xs.size() - 1));
  ArrayXd bids = Interpolate(bids_in, integrate_vals);
  ArrayXd win_probs =
      Interpolate(internal_bids_, other_highest_cdfs_[id], bids);
  ArrayXd exp_second_bid_given_win =
      Interpolate(internal_bids_, exp_value_funcs_[id], bids);
  ArrayXd profits = integrate_vals - exp_second_bid_given_win;
  ArrayXd revenues = exp_second_bid_given_win * win_probs;
  ArrayXd likelihoods =
      Interpolate(internal_values_, value_pdfs_[id], integrate_vals);
  return Areas(integrate_vals, revenues * likelihoods).sum();
}

float SecondPrice::GetValue(const Scatter& bids_in, int id) const {
  if (!pre_calculated_) {
    Precalculate();
  }

  ArrayXd integrate_vals =
      ArrayXd::LinSpaced((bids_in.xs.size() - 1) * 3, bids_in.xs(0),
                         bids_in.xs(bids_in.xs.size() - 1));
  ArrayXd bids = Interpolate(bids_in, integrate_vals);
  ArrayXd win_probs =
      Interpolate(internal_bids_, other_highest_cdfs_[id], bids);
  ArrayXd realized_value = integrate_vals * win_probs;
  ArrayXd likelihoods =
      Interpolate(internal_values_, value_pdfs_[id], integrate_vals);
  return Areas(integrate_vals, realized_value * likelihoods).sum();
}


void SecondPrice::Precalculate() const {
  for (int i = 0; i < n_players_; ++i) {
    std::vector<ArrayXd> other_cdfs;
    for (int j = 0; j < n_players_; ++j) {
      if (i != j) {
        other_cdfs.emplace_back(bid_cdfs_[j]);
      }
    }
    other_highest_cdfs_[i] =
        KthLowestOrderStatisticCDF(internal_bids_, other_cdfs, n_players_ - 1);
    exp_value_funcs_[i] =
        ExpectedValueFunctionCDF(internal_bids_, other_highest_cdfs_[i]);
  }
  pre_calculated_ = true;
}

}  // namespace auctions
