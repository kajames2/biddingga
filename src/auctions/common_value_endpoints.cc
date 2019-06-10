#include "auctions/common_value_endpoints.h"

#include <iostream>

#include <algorithm>

#include <boost/math/distributions/uniform.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include "numericaldists/distribution_ops.h"
#include "numericaldists/function_ops.h"
#include "numericaldists/order_statistic_ops.h"

using namespace boost::math;
using namespace numericaldists;
using namespace Eigen;

namespace auctions {

CommonValueEndpoints::CommonValueEndpoints(int n_bidders,
                                           Distribution value_dist,
                                           Distribution error_dist,
                                           int n_internal_samples,
                                           int value_integration_samples)
    : n_players_(n_bidders),
      pre_calculated_(false),
      others_bids_cdfs_(n_bidders),
      bid_funcs_(n_bidders) {
  utility_funcs_ = std::vector<std::function<double(double)>>(
      n_players_, [](double val) { return val; });
  prob_weight_funcs_ = std::vector<std::function<double(double)>>(
      n_players_, [](double prob) { return prob; });
  internal_bids_ = ArrayXd::LinSpaced(n_internal_samples,
                                      lower(value_dist) + lower(error_dist),
                                      upper(value_dist) + upper(error_dist));
  internal_values_ = ArrayXd::LinSpaced(value_integration_samples,
                                        lower(value_dist) + lower(error_dist),
                                        upper(value_dist) + upper(error_dist));
  internal_signals_ = ArrayXd::LinSpaced(n_internal_samples,
                                         lower(value_dist) + lower(error_dist),
                                         upper(value_dist) + upper(error_dist));
  ArrayXd internal_errors = ArrayXd::LinSpaced(
      n_internal_samples, lower(error_dist), upper(error_dist));
  ArrayXd error_pdf = internal_errors.unaryExpr(
      [&error_dist](double x) -> double { return pdf(error_dist, x); });
  ArrayXd value_pdf = internal_values_.unaryExpr(
      [&value_dist](double x) -> double { return pdf(value_dist, x); });
  ArrayXXd temp_joint = JointPDFIndependent(value_pdf, error_pdf);
  ArrayXXd value_mesh = GetXMesh(internal_values_, internal_errors.size());
  ArrayXXd error_mesh = GetYMesh(internal_errors, internal_values_.size());
  ArrayXXd signal_mesh = value_mesh + error_mesh;
  ArrayXXd value_cdf = TwoRandomVariableFunctionCDF(
      internal_values_, internal_errors, temp_joint, value_mesh, signal_mesh,
      internal_values_, internal_signals_);
  value_pdf_ = PDF2D(internal_values_, internal_signals_, value_cdf);
}

void CommonValueEndpoints::AcceptStrategy(Scatter bid_func, int id) {
  bid_funcs_[id] = std::move(bid_func);
  pre_calculated_ = false;
}

// rel_bid_func -- x-values are different precisions, y-values are the bid
// relative to mstar given the precision
float CommonValueEndpoints::GetFitness(const Scatter& bid_func, int id) const {
  if (!pre_calculated_) {
    Precalculate();
  }

  ArrayXd integration_signals =
      ArrayXd::LinSpaced((bid_func.xs.size() - 1) * 3, bid_func.xs(0),
                         bid_func.xs(bid_func.xs.size() - 1));
  ArrayXXd likelihoods =
      Interpolate2D(internal_values_, internal_signals_, value_pdf_,
                    internal_values_, integration_signals);

  ArrayXd bids = Interpolate(bid_func, integration_signals);
  ArrayXXd win_probs(integration_signals.size(), internal_values_.size());
  ArrayXXd profits(integration_signals.size(), internal_values_.size());
  for (int v = 0; v < internal_values_.size(); ++v) {
    win_probs.col(v) =
        Interpolate(internal_bids_, others_bids_cdfs_[id].col(v), bids);
    profits.col(v) = internal_values_[v] - bids;
  }

  ArrayXXd utils =
      profits.unaryExpr(utility_funcs_[id]) *
          win_probs.unaryExpr(prob_weight_funcs_[id]) +
      utility_funcs_[id](0) * (1 - win_probs).unaryExpr(prob_weight_funcs_[id]);
  return Areas2D(internal_values_, integration_signals, utils * likelihoods)
      .sum();
}

void CommonValueEndpoints::Precalculate() const {
  std::vector<ArrayXd> bid_sets;
  for (int i = 0; i < n_players_; ++i) {
    bid_sets.push_back(Interpolate(bid_funcs_[i], internal_signals_));
  }

  others_bids_cdfs_ = std::vector<ArrayXXd>(
      n_players_,
      ArrayXXd::Zero(internal_bids_.size(), internal_values_.size()));

#pragma omp parallel for
  for (int v = 0; v < internal_values_.size(); ++v) {
    std::vector<ArrayXd> cdfs(n_players_,
                              ArrayXd::Zero(internal_signals_.size()));
    for (int i = 0; i < n_players_; ++i) {
      if (value_pdf_.col(v).sum() > 0) {
        cdfs[i] = RandomVariableFunctionCDF(
            internal_signals_, value_pdf_.col(v), bid_sets[i], internal_bids_);
      }
    }
    for (int i = 0; i < n_players_; ++i) {
      std::vector<ArrayXd> other_cdfs;
      for (int j = 0; j < n_players_; ++j) {
        if (i != j) {
          other_cdfs.emplace_back(cdfs[j]);
        }
      }
      others_bids_cdfs_[i].col(v) =
          HighestOrderStatisticCDF(internal_bids_, other_cdfs);
    }
  }
  pre_calculated_ = true;
}

}  // namespace auctions
