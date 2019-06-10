#include "auctions/all_pay.h"

#include <algorithm>
#include <iostream>

#include "numericaldists/distribution.h"
#include "numericaldists/distribution_ops.h"
#include "numericaldists/function_ops.h"
#include "numericaldists/order_statistic_ops.h"

#include <omp.h>
#include <eigen3/Eigen/Core>
#include <iostream>

using namespace numericaldists;
using namespace Eigen;

namespace auctions {

AllPay::AllPay(const std::vector<Distribution>& value_dists,
               int n_internal_samples)
    : bid_cdfs_(value_dists.size()),
      n_players_(value_dists.size()),
      n_internal_samples_(n_internal_samples),
      max_bid_(upper(value_dists)) {
  utility_funcs_ = std::vector<std::function<double(double)>>(
      n_players_, [](double val) { return val; });
  prob_weight_funcs_ = std::vector<std::function<double(double)>>(
      n_players_, [](double prob) { return prob; });
  for (const auto& dist : value_dists) {
    ArrayXd values =
        ArrayXd::LinSpaced(n_internal_samples, lower(dist), upper(dist));
    ArrayXd pdfs =
        values.unaryExpr([&dist](double x) -> double { return pdf(dist, x); });
    value_pdfs_.push_back({values, pdfs});
  }
}

void AllPay::AcceptStrategy(const Scatter& bids, int id) {
  ArrayXd bid_range = ArrayXd::LinSpaced(
      n_internal_samples_, bids.ys.minCoeff(), bids.ys.maxCoeff());
  ArrayXd interp_bids = Interpolate(bids, value_pdfs_[id].xs);
  ArrayXd cdf = RandomVariableFunctionCDF(
      value_pdfs_[id].xs, value_pdfs_[id].ys, interp_bids, bid_range);
  bid_cdfs_[id] = {bid_range, cdf};
}


float AllPay::GetFitness(const Scatter& bids_in, int id) const {
  ArrayXd integrate_vals =
      ArrayXd::LinSpaced((bids_in.xs.size() - 1) * 3, bids_in.xs(0),
                         bids_in.xs(bids_in.xs.size() - 1));
  ArrayXd bids = Interpolate(bids_in, integrate_vals);
  ArrayXd profits = integrate_vals - bids;
  ArrayXd win_probs = ArrayXd::Ones(integrate_vals.size());
  for (int j = 0; j < n_players_; ++j) {
    if (j != id) {
      win_probs *= Interpolate(bid_cdfs_[j], bids);
    }
  }
  ArrayXd utils =
      profits.unaryExpr(utility_funcs_[id]) *
          win_probs.unaryExpr(prob_weight_funcs_[id]) +
      (-bids).unaryExpr(utility_funcs_[id]) * (1 - win_probs).unaryExpr(prob_weight_funcs_[id]);
  ArrayXd likelihoods = Interpolate(value_pdfs_[id], integrate_vals);
  return Areas(integrate_vals, utils * likelihoods).sum();
}


float AllPay::GetRevenue(const Scatter& bids_in, int id) const {
  ArrayXd integrate_vals =
      ArrayXd::LinSpaced((bids_in.xs.size() - 1) * 3, bids_in.xs(0),
                         bids_in.xs(bids_in.xs.size() - 1));
  ArrayXd bids = Interpolate(bids_in, integrate_vals);
  ArrayXd win_probs = ArrayXd::Ones(integrate_vals.size());
  for (int j = 0; j < n_players_; ++j) {
    if (j != id) {
      win_probs *= Interpolate(bid_cdfs_[j], bids);
    }
  }
  ArrayXd revenues = bids;
  ArrayXd likelihoods = Interpolate(value_pdfs_[id], integrate_vals);
  return Areas(integrate_vals, revenues * likelihoods).sum();
}


float AllPay::GetValue(const Scatter& bids_in, int id) const {
  ArrayXd integrate_vals =
      ArrayXd::LinSpaced((bids_in.xs.size() - 1) * 3, bids_in.xs(0),
                         bids_in.xs(bids_in.xs.size() - 1));
  ArrayXd bids = Interpolate(bids_in, integrate_vals);
  ArrayXd win_probs = ArrayXd::Ones(integrate_vals.size());
  for (int j = 0; j < n_players_; ++j) {
    if (j != id) {
      win_probs *= Interpolate(bid_cdfs_[j], bids);
    }
  }
  ArrayXd realized_values = integrate_vals * win_probs;
  ArrayXd likelihoods = Interpolate(value_pdfs_[id], integrate_vals);
  return Areas(integrate_vals, realized_values * likelihoods).sum();
}

}  // namespace auctions
