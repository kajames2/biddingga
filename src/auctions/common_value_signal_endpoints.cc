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
                                                  ArrayXd uncertainties) const {
  return Interpolate2D(relative_pdfs_[id].xs + internal_values_[value_index],
                       relative_pdfs_[id].ys, relative_pdfs_[id].zs, midpoints,
                       uncertainties);
}

CommonValueSignalEndpoints::CommonValueSignalEndpoints(
    uniform_distribution<> value_dist, uniform_distribution<> error_dist,
    std::vector<int> n_draws, int n_internal_samples,
    int value_integration_samples)
    : n_players_(n_draws.size()),
      n_internal_samples_(n_internal_samples),
      n_draws_(n_draws),
      pre_calculated_(false),
      bid_funcs_(n_draws.size()),
      one_draw_bids_(n_draws.size()),
      others_bids_cdfs_(n_draws.size()),
      relative_pdfs_(n_draws.size()),
      error_dist_(error_dist) {
  double min_error = error_dist.lower();
  double max_error = error_dist.upper();
  double min_value = value_dist.lower();
  double max_value = value_dist.upper();
  double min_midpoint = min_value + min_error;
  double max_midpoint = max_value + max_error;
  double error_range = max_error - min_error;
  internal_bids_ =
      ArrayXd::LinSpaced(n_internal_samples_, min_midpoint, max_midpoint);
  internal_values_ =
      ArrayXd::LinSpaced(value_integration_samples, min_midpoint, max_midpoint);
  internal_midpoints_ =
      ArrayXd::LinSpaced(n_internal_samples_, min_midpoint, max_midpoint);
  internal_uncertainties_ =
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
  ArrayXd rel_midpoints = ArrayXd::LinSpaced(
      n_internal_samples_, min_error - buffer, max_error + buffer);
  ArrayXd buffered_unc =
      ArrayXd::LinSpaced(n_internal_samples_, 0, error_range + 2 * buffer);
  ArrayXd rel_midpoint_cdf = rel_midpoints.unaryExpr(
      [&error_dist](double x) -> double { return cdf(error_dist, x); });
  ArrayXXd min_midpoint_mesh = GetXMesh(rel_midpoints, n_internal_samples_);
  ArrayXXd max_midpoint_mesh = GetYMesh(rel_midpoints, n_internal_samples_);
  ArrayXXd midpoint_mesh = (min_midpoint_mesh + max_midpoint_mesh) / 2;
  ArrayXXd unc_mesh =
      error_range - (max_midpoint_mesh - min_midpoint_mesh).abs();
  for (int i = 0; i < n_players_; ++i) {
    int draw = n_draws_[i];
    if (draw > 1) {
      ArrayXXd joint = LowestHighestJointOrderStatisticPDF(
          rel_midpoints, rel_midpoint_cdf, draw);
      ArrayXXd joint_m_r_cdf = TwoRandomVariableFunctionCDF(
          rel_midpoints, rel_midpoints, joint, midpoint_mesh, unc_mesh,
          rel_midpoints, buffered_unc);
      relative_pdfs_[i] = {rel_midpoints, buffered_unc,
                           PDF2D(rel_midpoints, buffered_unc, joint_m_r_cdf)};
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

// bid_func -- x-values are different uncertainties, y-values are the bid
// relative to midpoint given the uncertainty
float CommonValueSignalEndpoints::GetFitness(const Grid& bid_func,
                                             int id) const {
  if (!pre_calculated_) {
    Precalculate();
  }
  ArrayXd integration_midpoints =
      ArrayXd::LinSpaced((bid_func.xs.size() - 1) * 3, bid_func.xs(0),
                         bid_func.xs(bid_func.xs.size() - 1));
  ArrayXd integration_uncertainties =
      ArrayXd::LinSpaced((bid_func.ys.size() - 1) * 3, bid_func.ys(0),
                         bid_func.ys(bid_func.ys.size() - 1));
  double min_integration_value = integration_midpoints(0) + error_dist_.lower();
  double max_integration_value =
      integration_midpoints(integration_midpoints.size() - 1) +
      error_dist_.upper();
  ArrayXd value_likelihoods = value_pdf_;
  ArrayXXd bids =
      Interpolate2D(bid_func, integration_midpoints, integration_uncertainties);
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
    ArrayXXd win_probs = ArrayXXd::Zero(integration_midpoints.size(),
                                        integration_uncertainties.size());
    for (int i = 0; i < bids.rows(); ++i) {
      win_probs.row(i) = Interpolate(internal_bids_,
                                     others_bids_cdfs_[id].col(v), bids.row(i));
    }
    ArrayXXd signal_likelihoods =
        GetSignalPDF(id, v, integration_midpoints, integration_uncertainties);
    ArrayXXd likelihoods = value_likelihoods[v] * signal_likelihoods;
    ArrayXXd profits = internal_values_[v] - bids;
    ArrayXXd utils = profits.unaryExpr(utility_funcs_[id]) *
                         win_probs.unaryExpr(prob_weight_funcs_[id]) +
                     utility_funcs_[id](0) *
                         (1 - win_probs).unaryExpr(prob_weight_funcs_[id]);
    sub_areas(v) = Areas2D(integration_midpoints, integration_uncertainties,
                           utils * likelihoods)
                       .sum();
  }
  return Areas(internal_values_, sub_areas).sum();
}

float CommonValueSignalEndpoints::GetFitness(const Scatter& bid_func,
                                             int id) const {
  if (!pre_calculated_) {
    Precalculate();
  }
  ArrayXd integration_midpoints =
      ArrayXd::LinSpaced((bid_func.xs.size() - 1) * 3, bid_func.xs(0),
                         bid_func.xs(bid_func.xs.size() - 1));
  double min_integration_value = integration_midpoints(0) + error_dist_.lower();
  double max_integration_value =
      integration_midpoints(integration_midpoints.size() - 1) +
      error_dist_.upper();
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
  std::vector<ArrayXXd> bid_sets(n_players_);
  std::vector<ArrayXd> one_draw_bid_sets(n_players_);
  for (int i = 0; i < n_players_; ++i) {
    if (n_draws_[i] > 1) {
      bid_sets[i] = Interpolate2D(bid_funcs_[i], internal_midpoints_,
                                  internal_uncertainties_);
    } else {
      one_draw_bid_sets[i] =
          Interpolate(one_draw_bids_[i], internal_midpoints_);
    }
  }
  others_bids_cdfs_ = std::vector<ArrayXXd>(
      n_players_, ArrayXXd::Zero(n_internal_samples_, n_internal_samples_));
#pragma omp parallel for
  for (int v = 0; v < n_internal_samples_; ++v) {
    std::vector<ArrayXd> cdfs(n_players_, ArrayXd::Zero(n_internal_samples_));
    if (value_pdf_(v) > 0) {
      for (int i = 0; i < n_players_; ++i) {
        if (n_draws_[i] > 1) {
          ArrayXXd signal_likelihoods =
              GetSignalPDF(i, v, internal_midpoints_, internal_uncertainties_);
          cdfs[i] = RandomVariableFunctionCDF(
              internal_midpoints_, internal_uncertainties_, signal_likelihoods,
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
  pre_calculated_ = true;
}

}  // namespace auctions
