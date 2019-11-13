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

ArrayXXd CommonValueSignalEndpoints::GetSignalPDF(int id, int value_index,
                                                  ArrayXd midpoints,
                                                  ArrayXd precisions) const {
  auto pdf = ArrayXXd::Zero(midpoints.size(), precisions.size());
  Interpolate2D(relative_pdfs_[id].xs + internal_values_[value_index],
                relative_pdfs_[id].ys, relative_pdfs_[id].zs, midpoints,
                precisions);
  return pdf;
}

CommonValueSignalEndpoints::CommonValueSignalEndpoints(
    Distribution value_dist, Distribution error_dist, std::vector<int> n_draws,
    int n_internal_samples, int value_integration_samples)
    : n_players_(n_draws.size()),
      n_internal_samples_(n_internal_samples),
      n_draws_(n_draws),
      pre_calculated_(false),
      bid_funcs_(n_draws.size()),
      one_draw_bids_(n_draws.size()),
      others_bids_cdfs_(n_draws.size()),
      relative_pdfs_(n_draws.size()),
      error_dist_(error_dist) {
  double min_error = lower(error_dist);
  double max_error = upper(error_dist);
  double min_value = lower(value_dist);
  double max_value = upper(value_dist);
  double min_midpoint = min_value + min_error;
  double max_midpoint = max_value + max_error;
  double error_range = max_error - min_error;
  internal_bids_ =
      ArrayXd::LinSpaced(n_internal_samples_, min_midpoint, max_midpoint);
  internal_values_ =
      ArrayXd::LinSpaced(value_integration_samples, min_midpoint, max_midpoint);
  internal_midpoints_ =
      ArrayXd::LinSpaced(n_internal_samples_, min_midpoint, max_midpoint);
  internal_precisions_ =
      ArrayXd::LinSpaced(n_internal_samples_, 0, error_range);

  ArrayXd internal_errors =
      ArrayXd::LinSpaced(n_internal_samples_, min_error, max_error);
  value_pdf_ = internal_values_.unaryExpr(
      [&value_dist](double x) -> double { return pdf(value_dist, x); });
  utility_funcs_ = std::vector<std::function<double(double)>>(
      n_players_, [](double val) { return val; });
  prob_weight_funcs_ = std::vector<std::function<double(double)>>(
      n_players_, [](double prob) { return prob; });

  double buffer = 0.05 * error_range;
  auto rel_midpoints = ArrayXd::LinSpaced(
      n_internal_samples_, min_error - buffer, max_error + buffer);
  auto buffered_precs =
      ArrayXd::LinSpaced(n_internal_samples_, 0, error_range + 2 * buffer);
  ArrayXd rel_midpoint_cdf = rel_midpoints.unaryExpr(
      [&error_dist](double x) -> double { return cdf(error_dist, x); });
  ArrayXXd min_midpoint_mesh = GetXMesh(rel_midpoints, n_internal_samples_);
  ArrayXXd max_midpoint_mesh = GetYMesh(rel_midpoints, n_internal_samples_);
  ArrayXXd midpoint_mesh = (min_midpoint_mesh + max_midpoint_mesh) / 2;
  ArrayXXd prec_mesh = (max_midpoint_mesh - min_midpoint_mesh).abs();
  for (int i = 0; i < n_players_; ++i) {
    int draw = n_draws_[i];
    if (draw > 1) {
      ArrayXXd joint = LowestHighestJointOrderStatisticPDF(
          rel_midpoints, rel_midpoint_cdf, draw);
      ArrayXXd joint_m_r_cdf = TwoRandomVariableFunctionCDF(
          rel_midpoints, rel_midpoints, joint, midpoint_mesh, prec_mesh,
          internal_midpoints_, buffered_precs);
      relative_pdfs_[i] = {rel_midpoints, buffered_precs,
                           PDF2D(rel_midpoints, buffered_precs, joint_m_r_cdf)};
    }
  }
}

void CommonValueSignalEndpoints::AcceptStrategy(Grid bid_func, int id) {
  bid_funcs_[id] = std::move(bid_func);
  pre_calculated_ = false;
}

void CommonValueSignalEndpoints::AcceptStrategy(Scatter bid_func, int id) {
  one_draw_bids_[id] = std::move(bid_func);
  pre_calculated_ = false;
}

// bid_func -- x-values are different precisions, y-values are the bid
// relative to midpoint given the precision
float CommonValueSignalEndpoints::GetFitness(const Grid& bid_func,
                                             int id) const {
  if (!pre_calculated_) {
    Precalculate();
  }
  ArrayXd integration_midpoints =
      ArrayXd::LinSpaced((bid_func.xs.size() - 1) * 3, bid_func.xs(0),
                         bid_func.xs(bid_func.xs.size() - 1));
  ArrayXd integration_precisions =
      ArrayXd::LinSpaced((bid_func.ys.size() - 1) * 3, bid_func.ys(0),
                         bid_func.ys(bid_func.ys.size() - 1));
  double min_integration_value = integration_midpoints(0) + lower(error_dist_);
  double max_integration_value =
      integration_midpoints(integration_midpoints.size() - 1) +
      upper(error_dist_);
  ArrayXd value_likelihoods = value_pdf_;
  ArrayXXd bids =
      Interpolate2D(bid_func, integration_midpoints, integration_precisions);
  ArrayXd sub_areas = ArrayXd::Zero(internal_values_.size());
  int v = 0;
  while (v < internal_values_.size() - 1 &&
         internal_values_(v + 1) < min_integration_value) {
    ++v;
  }
  for (; v < internal_values_.size(); ++v) {
    if (v > 0 && internal_values_(v - 1) > max_integration_value) {
      break;
    }
    auto win_probs = ArrayXXd::Zero(integration_midpoints.size(),
                                    integration_precisions.size());
    for (int i = 0; i < bids.rows(); ++i) {
      win_probs.row(i) = Interpolate(internal_bids_,
                                     others_bids_cdfs_[id].col(v), bids.row(i));
    }
    ArrayXXd signal_likelihoods =
        GetSignalPDF(id, v, integration_midpoints, integration_precisions);
    ArrayXXd likelihoods = value_likelihoods[v] * signal_likelihoods;
    ArrayXXd profits = internal_values_[v] - bids;
    ArrayXXd utils = profits.unaryExpr(utility_funcs_[id]) *
                         win_probs.unaryExpr(prob_weight_funcs_[id]) +
                     utility_funcs_[id](0) *
                         (1 - win_probs).unaryExpr(prob_weight_funcs_[id]);
    sub_areas(v) = Areas2D(integration_midpoints, integration_precisions,
                           utils * likelihoods[v])
                       .sum();
  }
  return Areas(internal_values_, sub_areas).sum();
}  // namespace auctions

float CommonValueSignalEndpoints::GetFitness(const Scatter& bid_func,
                                             int id) const {
  if (!pre_calculated_) {
    Precalculate();
  }
  ArrayXd integration_midpoints =
      ArrayXd::LinSpaced((bid_func.xs.size() - 1) * 3, bid_func.xs(0),
                         bid_func.xs(bid_func.xs.size() - 1));
  double min_integration_value = integration_midpoints(0) + lower(error_dist_);
  double max_integration_value =
      integration_midpoints(integration_midpoints.size() - 1) +
      upper(error_dist_);
  ArrayXd value_likelihoods = value_pdf_;
  ArrayXXd win_probs =
      ArrayXXd::Zero(integration_midpoints.size(), internal_values_.size());
  ArrayXXd profits =
      ArrayXXd::Zero(integration_midpoints.size(), internal_values_.size());
  ArrayXXd likelihoods =
      ArrayXXd::Zero(integration_midpoints.size(), internal_values_.size());
  ArrayXd bids = Interpolate(bid_func, integration_midpoints);
  int v = 0;
  while (v < internal_values_.size() - 1 &&
         internal_values_(v + 1) < min_integration_value) {
    ++v;
  }
  for (; v < internal_values_.size(); ++v) {
    if (v > 0 && internal_values_(v - 1) > max_integration_value) {
      break;
    }
    ArrayXd error_likelihoods = (integration_midpoints - internal_values_[v])
                                    .unaryExpr([this](double x) -> double {
                                      return pdf(error_dist_, x);
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
  return Areas2D(internal_values_, integration_midpoints, utils * likelihoods)
      .sum();
}

void CommonValueSignalEndpoints::Precalculate() const {
  std::vector<ArrayXXd> bid_sets;
  std::vector<ArrayXd> one_draw_bid_sets;
  for (const auto& bid : bid_funcs_) {
    bid_sets.push_back(
        Interpolate2D(bid, internal_midpoints_, internal_precisions_));
  }
  for (const auto& bid : one_draw_bids_) {
    one_draw_bid_sets.push_back(Interpolate(bid, internal_midpoints_));
  }
  others_bids_cdfs_ = std::vector<ArrayXXd>(
      n_players_, ArrayXXd::Zero(n_internal_samples_, n_internal_samples_));
#pragma omp parallel for
  for (int v = 0; v < n_internal_samples_; ++v) {
    std::vector<ArrayXd> cdfs(n_players_, ArrayXd::Zero(n_internal_samples_));
    if (value_pdf_(v) > 0) {
      for (int i = 0; i < n_players_; ++i) {
        if (n_draws_[i] > 1) {
          auto signal_likelihoods =
              GetSignalPDF(i, v, internal_midpoints_, internal_precisions_);
          cdfs[i] = RandomVariableFunctionCDF(
              internal_midpoints_, internal_precisions_, signal_likelihoods,
              bid_sets[i], internal_bids_);
        } else {
          ArrayXd error_likelihoods =
              (internal_midpoints_ - v).unaryExpr([this](double x) -> double {
                return pdf(error_dist_, x);
              });
          cdfs[i] =
              RandomVariableFunctionCDF(internal_midpoints_, error_likelihoods,
                                        one_draw_bid_sets[i], internal_bids_);
        }
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
}

}  // namespace auctions
