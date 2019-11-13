#include "auctions/common_value_signal_endpoints.h"

#include <boost/math/distributions/uniform.hpp>
#include <cstdlib>
#include <vector>
#include "numericaldists/distribution_ops.h"
#include "numericaldists/function_ops.h"
#include "numericaldists/order_statistic_ops.h"

using namespace boost::math;
using namespace numericaldists;
using namespace Eigen;

namespace auctions {

ArrayXXd GetSignalPDF(int value_index) {
  pdf = ArrayXXd::Zero(internal_precisions_, internal_midpoints_);
  pdf.block(0, value_index, internal_precisions_.size(), error_size_) =
      joint_pdf;
  return pdf;
}

CommonValueSignalEndpoints::CommonValueSignalEndpoints(
    Distribution value_dist, Distribution error_dist, std::vector<int> n_draws,
    int n_internal_samples, int value_integration_samples)
    : n_players_(n_draws.size()),
      pre_calculated_(false),
      bid_funcs_(n_draws.size()),
      n_draws_(n_draws),
      one_draw_bids_(n_draws.size()),
      others_bids_cdfs_(n_draws.size()),
      value_dists_(n_draws.size()),
      n_internal_samples_(n_internal_samples),
      mstar_integration_samples_(mstar_integration_samples) {
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
  utility_funcs_ = std::vector<std::function<double(double)>>(
      n_players_, [](double val) { return val; });
  prob_weight_funcs_ = std::vector<std::function<double(double)>>(
      n_players_, [](double prob) { return prob; });
  internal_bids_ = ArrayXd::LinSpaced(n_internal_samples_, -epsilon +
                                      : id_int.min, epsilon + bid_int.max);
  internal_signals_ =
      ArrayXd::LinSpaced(n_internal_samples_, -epsilon, epsilon);
  internal_mstars_ = ArrayXd::LinSpaced(n_internal_samples_, -epsilon, epsilon);
  internal_precs_ = ArrayXd::LinSpaced(n_internal_samples_, 0, 2 * epsilon);
  one_draw_pdf_ = ArrayXd::Ones(n_internal_samples, 1) / (2 * epsilon);
  ArrayXd signal_cdf = ArrayXd::LinSpaced(n_internal_samples_, 0, 1);
  ArrayXXd signalXMesh = GetXMesh(internal_signals_, internal_signals_.size());
  ArrayXXd signalYMesh = GetYMesh(internal_signals_, internal_signals_.size());
  ArrayXXd mstarMesh = (signalXMesh + signalYMesh) / 2;
  ArrayXXd precMesh = (signalYMesh - signalXMesh).abs();
  for (int i = 0; i < n_players_; ++i) {
    int draw = n_draws_[i];
    if (draw > 1) {
      ArrayXXd joint = LowestHighestJointOrderStatisticPDF(internal_signals_,
                                                           signal_cdf, draw);
      ArrayXXd joint_m_r_cdf = TwoRandomVariableFunctionCDF(
          internal_signals_, internal_signals_, joint, mstarMesh, precMesh,
          internal_mstars_, internal_precs_);
      value_dists_[i] = PDF2D(internal_mstars_, internal_precs_, joint_m_r_cdf);
    }
  }
}

void CommonValueSignalEndpoints::AcceptStrategy(Grid bid_func, int id) {
  bid_funcs_[id] = std::move(bid_func);
  pre_calculated_ = false;
}

void CommonValueSignalEndpoints::AcceptStrategy(Scatter bid_func int id) {
  one_draw_bids_[id] = bid_func;
  pre_calculated_ = false;
}

// bid_func -- x-values are different precisions, y-values are the bid
// relative to mstar given the precision
float CommonValueSignalEndpoints::GetFitness(const Grid& bid_func,
                                             int id) const {
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
  for (int v = 0; v < internal_values_.size(); ++v) {
    ArrayXd bids = Interpolate(bid_func, integration_signals);
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
  if (!pre_calculated_) {
    Precalculate();
  }

  ArrayXd integrate_mstars =
      ArrayXd::LinSpaced(mstar_integration_samples_, internal_mstars_(0),
                         internal_mstars_(internal_mstars_.size() - 1));
  ArrayXd integrate_precs =
      ArrayXd::LinSpaced((bid_func.xs.size() - 1) * 3, bid_func.xs.minCoeff(),
                         bid_func.xs.maxCoeff());

  ArrayXXd mstar_mesh = GetXMesh(integrate_mstars, integrate_precs.size());
  ArrayXXd prec_mesh = GetYMesh(integrate_precs, integrate_mstars.size());
  ArrayXd bids = Interpolate(bid_func, integrate_precs);
  ArrayXXd bids_mesh = GetYMesh(bids, integrate_mstars.size());
  ArrayXXd bids = mstar_mesh + bids_mesh;
  ArrayXXd win_probs(bids_mesh.rows(), bids_mesh.cols());
  for (int i = 0; i < bids_mesh.rows(); ++i) {
    win_probs.row(i) =
        Interpolate(internal_bids_, others_bids_cdfs_[id], bids.row(i));
  }
  // True value is the reference point, so it is 0.
  ArrayXXd profits = 0 - bids;
  ArrayXXd utils =
      profits.unaryExpr(utility_funcs_[id]) *
          win_probs.unaryExpr(prob_weight_funcs_[id]) +
      utility_funcs_[id](0) * (1 - win_probs).unaryExpr(prob_weight_funcs_[id]);
  ArrayXXd likelihoods =
      Interpolate2D(internal_mstars_, internal_precs_, value_dists_[id],
                    integrate_mstars, integrate_precs);
  return Areas2D(integrate_mstars, integrate_precs, utils * likelihoods).sum();
}

float CommonValueSignalEndpoints::GetFitness(const float bid, int id) const {
  if (!pre_calculated_) {
    Precalculate();
  }

  ArrayXd integrate_signals =
      ArrayXd::LinSpaced(mstar_integration_samples_, internal_signals_(0),
                         internal_signals_(internal_signals_.size() - 1));
  ArrayXd bids = integrate_signals + bid;
  ArrayXd win_probs = Interpolate(internal_bids_, others_bids_cdfs_[id], bids);
  ArrayXd profits = 0 - bids;
  ArrayXd utils =
      profits.unaryExpr(utility_funcs_[id]) *
          win_probs.unaryExpr(prob_weight_funcs_[id]) +
      utility_funcs_[id](0) * (1 - win_probs).unaryExpr(prob_weight_funcs_[id]);
  ArrayXd likelihoods =
      Interpolate(internal_signals_, one_draw_pdf_, integrate_signals);
  return Areas(integrate_signals, utils * likelihoods).sum();
}

void CommonValueSignalEndpoints::Precalculate() const {
  std::vector<ArrayXXd> bid_sets;
  for (int i = 0; i < n_players_; ++i) {
    bid_sets.push_back(Interpolate2D(bid_funcs_[i], internal_midpoints_,
                                     internal_precisions_));
  }
  others_bids_cdfs_ = std::vector<ArrayXXd>(
      n_players_,
      ArrayXXd::Zero(internal_bids_.size(), internal_values_.size()));
#pragma omp parallel for
  for (int v = 0; v < internal_values_.size(); ++v) {
    signal_likelihoods = GetSignalPDF(v);
    std::vector<ArrayXd> cdfs(n_players_,
                              ArrayXd::Zero(internal_signals_.size()));
    if (value_pdf_(v) > 0) {
      for (int i = 0; i < n_players_; ++i) {
        cdfs[i] = RandomVariableFunctionCDF(
            internal_midpoints_, internal_precisions_, signal_likelihoods,
            bid_sets[i], internal_bids_);
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
