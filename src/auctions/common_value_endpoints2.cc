#include "auctions/common_value_endpoints2.h"

#include <boost/math/distributions/uniform.hpp>
#include <vector>
#include "numericaldists/distribution_ops.h"
#include "numericaldists/function_ops.h"
#include "numericaldists/order_statistic_ops.h"

using namespace boost::math;
using namespace numericaldists;
using namespace Eigen;

namespace auctions {

CommonValueEndpoints2::CommonValueEndpoints2(int n_bidders,
                                             Distribution value_dist,
                                             Distribution error_dist,
                                             int n_internal_samples,
                                             int value_integration_samples)
    : n_players_(n_bidders),
      pre_calculated_(false),
      others_bids_cdfs_(n_bidders),
      bid_funcs_(n_bidders),
      error_dist_(error_dist) {
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
  value_pdf_ = internal_values_.unaryExpr(
      [&value_dist](double x) -> double { return pdf(value_dist, x); });
}

void CommonValueEndpoints2::AcceptStrategy(Scatter bid_func, int id) {
  bid_funcs_[id] = std::move(bid_func);
  pre_calculated_ = false;
}

float CommonValueEndpoints2::GetFitness(const Scatter& bid_func, int id) const {
  if (!pre_calculated_) {
    Precalculate();
  }
  ArrayXd integration_signals =
      ArrayXd::LinSpaced((bid_func.xs.size() - 1) * 3, bid_func.xs(0),
                         bid_func.xs(bid_func.xs.size() - 1));
  ArrayXd value_likelihoods = value_pdf_;
  ArrayXXd win_probs(integration_signals.size(), internal_values_.size());
  ArrayXXd profits(integration_signals.size(), internal_values_.size());
  ArrayXXd likelihoods(integration_signals.size(), internal_values_.size());
  Distribution error_dist = error_dist_;
  ArrayXd bids = Interpolate(bid_func, integration_signals);
  for (int v = 0; v < internal_values_.size(); ++v) {
    ArrayXd error_likelihoods =
        (integration_signals - internal_values_[v])
            .unaryExpr([&error_dist](double x) -> double {
              return pdf(error_dist, x);
            });
    likelihoods.col(v) = value_likelihoods[v] * error_likelihoods;
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

void CommonValueEndpoints2::Precalculate() const {
  std::vector<ArrayXd> bid_sets;
  for (int i = 0; i < n_players_; ++i) {
    bid_sets.push_back(Interpolate(bid_funcs_[i], internal_signals_));
  }

  others_bids_cdfs_ = std::vector<ArrayXXd>(
      n_players_,
      ArrayXXd::Zero(internal_bids_.size(), internal_values_.size()));
  Distribution error_dist = error_dist_;
#pragma omp parallel for
  for (int v = 0; v < internal_values_.size(); ++v) {
    ArrayXd error_likelihoods =
        (internal_signals_ - internal_values_[v])
            .unaryExpr([&error_dist](double x) -> double {
              return pdf(error_dist, x);
            });
    std::vector<ArrayXd> cdfs(n_players_,
                              ArrayXd::Zero(internal_signals_.size()));
    if (value_pdf_(v) > 0) {
      for (int i = 0; i < n_players_; ++i) {
        cdfs[i] = RandomVariableFunctionCDF(
            internal_signals_, error_likelihoods, bid_sets[i], internal_bids_);
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
