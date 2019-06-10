#include "auctions/first_price_reverse.h"

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

FirstPriceReverse::FirstPriceReverse(
    const std::vector<Distribution>& cost_dists, int n_internal_samples)
    : bid_cdfs_(cost_dists.size()),
      n_players_(cost_dists.size()),
      n_internal_samples_(n_internal_samples) {
  utility_funcs_ = std::vector<std::function<double(double)>>(
      n_players_, [](double val) { return val; });
  prob_weight_funcs_ = std::vector<std::function<double(double)>>(
      n_players_, [](double prob) { return prob; });
  for (const auto& dist : cost_dists) {
    ArrayXd values =
        ArrayXd::LinSpaced(n_internal_samples, lower(dist), upper(dist));
    ArrayXd pdfs =
        values.unaryExpr([&dist](double x) -> double { return pdf(dist, x); });
    cost_pdfs_.push_back({values, pdfs});
  }
}

void FirstPriceReverse::AcceptStrategy(const Scatter& bids, int id) {
  ArrayXd bid_range = ArrayXd::LinSpaced(
      n_internal_samples_, bids.ys.minCoeff(), bids.ys.maxCoeff());
  ArrayXd interp_bids = Interpolate(bids, cost_pdfs_[id].xs);
  ArrayXd cdf = RandomVariableFunctionCDF(
      cost_pdfs_[id].xs, cost_pdfs_[id].ys, interp_bids, bid_range);
  bid_cdfs_[id] = {bid_range, cdf};
}

float FirstPriceReverse::GetFitness(const Scatter& bids_in, int id) const {
  ArrayXd integrate_costs = ArrayXd::LinSpaced(
      (bids_in.xs.size() - 1) * 3, bids_in.xs(0), bids_in.xs(bids_in.xs.size() - 1));
  ArrayXd bids = Interpolate(bids_in, integrate_costs);
  ArrayXd profits = bids - integrate_costs;
  ArrayXd win_probs = ArrayXd::Ones(integrate_costs.size());
  for (int j = 0; j < n_players_; ++j) {
    if (j != id) {
      win_probs *= 1 - Interpolate(bid_cdfs_[j], bids);
    }
  }
  ArrayXd utils =
      profits.unaryExpr(utility_funcs_[id]) *
          win_probs.unaryExpr(prob_weight_funcs_[id]) +
      utility_funcs_[id](0) * (1 - win_probs).unaryExpr(prob_weight_funcs_[id]);
  ArrayXd likelihoods = Interpolate(cost_pdfs_[id], integrate_costs);
  return Areas(integrate_costs, utils * likelihoods).sum();
}

}  // namespace auctions
