#include "auctions/common_value_signal.h"

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

CommonValueSignal::CommonValueSignal(std::vector<int> n_draws, float epsilon,
                                     Interval rel_bid_int,
                                     int n_internal_samples,
                                     int mstar_integration_samples)
    : n_players_(n_draws.size()),
      pre_calculated_(false),
      n_draws_(n_draws),
      rel_bid_funcs_(n_draws.size()),
      one_draw_rel_bids_(n_draws.size()),
      others_bids_cdfs_(n_draws.size()),
      value_dists_(n_draws.size()),
      n_internal_samples_(n_internal_samples),
      mstar_integration_samples_(mstar_integration_samples) {
  utility_funcs_ = std::vector<std::function<double(double)>>(
      n_players_, [](double val) { return val; });
  prob_weight_funcs_ = std::vector<std::function<double(double)>>(
      n_players_, [](double prob) { return prob; });
  internal_bids_ =
      ArrayXd::LinSpaced(n_internal_samples_, -epsilon + rel_bid_int.min,
                         epsilon + rel_bid_int.max);
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

void CommonValueSignal::AcceptStrategy(Scatter rel_bid_func, int id) {
  rel_bid_funcs_[id] = std::move(rel_bid_func);
  pre_calculated_ = false;
}

void CommonValueSignal::AcceptStrategy(float rel_bid_value, int id) {
  one_draw_rel_bids_[id] = rel_bid_value;
  pre_calculated_ = false;
}

// rel_bid_func -- x-values are different precisions, y-values are the bid
// relative to mstar given the precision
float CommonValueSignal::GetFitness(const Scatter& rel_bid_func, int id) const {
  if (!pre_calculated_) {
    Precalculate();
  }

  ArrayXd integrate_mstars =
      ArrayXd::LinSpaced(mstar_integration_samples_, internal_mstars_(0),
                         internal_mstars_(internal_mstars_.size() - 1));
  ArrayXd integrate_precs = ArrayXd::LinSpaced((rel_bid_func.xs.size() - 1) * 3,
                                               rel_bid_func.xs.minCoeff(),
                                               rel_bid_func.xs.maxCoeff());

  ArrayXXd mstar_mesh = GetXMesh(integrate_mstars, integrate_precs.size());
  ArrayXXd prec_mesh = GetYMesh(integrate_precs, integrate_mstars.size());
  ArrayXd rel_bids = Interpolate(rel_bid_func, integrate_precs);
  ArrayXXd rel_bids_mesh = GetYMesh(rel_bids, integrate_mstars.size());
  ArrayXXd bids = mstar_mesh + rel_bids_mesh;
  ArrayXXd win_probs(rel_bids_mesh.rows(), rel_bids_mesh.cols());
  for (int i = 0; i < rel_bids_mesh.rows(); ++i) {
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

float CommonValueSignal::GetFitness(const float rel_bid, int id) const {
  if (!pre_calculated_) {
    Precalculate();
  }

  ArrayXd integrate_signals =
      ArrayXd::LinSpaced(mstar_integration_samples_, internal_signals_(0),
                         internal_signals_(internal_signals_.size() - 1));
  ArrayXd bids = integrate_signals + rel_bid;
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

void CommonValueSignal::Precalculate() const {
  std::vector<ArrayXd> cdfs(n_players_);

  #pragma omp parallel for
  for (int i = 0; i < n_players_; ++i) {
    if (n_draws_[i] > 1) {
      ArrayXd internal_rel_bids =
          Interpolate(rel_bid_funcs_[i], internal_precs_);
      ArrayXXd rel_bid_mesh =
          GetYMesh(internal_rel_bids, internal_mstars_.size());
      ArrayXXd mstar_mesh =
          GetXMesh(internal_mstars_, internal_rel_bids.size());
      ArrayXXd bid_mesh = mstar_mesh + rel_bid_mesh;
      cdfs[i] =
          RandomVariableFunctionCDF(internal_mstars_, internal_precs_,
                                    value_dists_[i], bid_mesh, internal_bids_);
    } else {
      ArrayXd bids = internal_mstars_ + one_draw_rel_bids_[i];
      cdfs[i] = RandomVariableFunctionCDF(internal_mstars_, one_draw_pdf_, bids,
                                          internal_bids_);
    }
  }

  for (int i = 0; i < n_players_; ++i) {
    std::vector<ArrayXd> other_cdfs;
    for (int j = 0; j < n_players_; ++j) {
      if (i != j) {
        other_cdfs.emplace_back(cdfs[j]);
      }
    }
    others_bids_cdfs_[i] =
        KthLowestOrderStatisticCDF(internal_bids_, other_cdfs, n_players_ - 1);
  }
  pre_calculated_ = true;
}

}  // namespace auctions
