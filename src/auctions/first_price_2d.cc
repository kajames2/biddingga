#include "auctions/first_price_2d.h"

#include <algorithm>
#include <iostream>

#include "numericaldists/distribution.h"
#include "numericaldists/distribution_ops.h"
#include "numericaldists/function_ops.h"
#include "numericaldists/grid_multi.h"
#include "numericaldists/order_statistic_ops.h"

#include <omp.h>
#include <eigen3/Eigen/Core>
#include <iostream>

using namespace numericaldists;
using namespace Eigen;

namespace auctions {

void FirstPrice2D::Initialize(const std::vector<Distribution>& valuex_dists,
                              const std::vector<Distribution>& valuey_dists) {
  prob_weight_funcs_ = std::vector<std::function<double(double)>>(
      n_players_, [](double prob) { return prob; });
  for (int i = 0; i < n_players_; ++i) {
    auto distx = valuex_dists[i];
    auto disty = valuey_dists[i];
    ArrayXd valuesx =
        ArrayXd::LinSpaced(n_internal_samples_, lower(distx), upper(distx));
    ArrayXd valuesy =
        ArrayXd::LinSpaced(n_internal_samples_, lower(disty), upper(disty));
    ArrayXd xpdfs = valuesx.unaryExpr(
        [&distx](double x) -> double { return pdf(distx, x); });
    ArrayXd ypdfs = valuesy.unaryExpr(
        [&disty](double y) -> double { return pdf(disty, y); });
    value_pdfs_.push_back(
        {valuesx, valuesy, JointPDFIndependent(xpdfs, ypdfs)});
  }
}

FirstPrice2D::FirstPrice2D(const std::vector<Distribution>& valuex_dists,
                           const std::vector<Distribution>& valuey_dists,
                           int n_internal_samples)
    : bid_cdfs_(valuex_dists.size()),
      n_players_(valuex_dists.size()),
      bid_margx_cdfs_(valuex_dists.size()),
      bid_margy_cdfs_(valuex_dists.size()),
      n_internal_samples_(n_internal_samples) {
  std::function<double(double, double)> util = [](double x, double y) {
    return x + y;
  };
  utility_funcs_ =
      std::vector<std::function<double(double, double)>>(n_players_, util);
  Initialize(valuex_dists, valuey_dists);
}

FirstPrice2D::FirstPrice2D(
    const std::vector<Distribution>& valuex_dists,
    const std::vector<Distribution>& valuey_dists,
    std::vector<std::function<double(double, double)>> utils,
    int n_internal_samples)
    : bid_cdfs_(valuex_dists.size()),
      n_players_(valuex_dists.size()),
      bid_margx_cdfs_(valuex_dists.size()),
      bid_margy_cdfs_(valuex_dists.size()),
      n_internal_samples_(n_internal_samples),
      utility_funcs_(utils) {
  Initialize(valuex_dists, valuey_dists);
}

void FirstPrice2D::AcceptStrategy(const GridMulti& bids, int id) {
  ArrayXd bidx_range =
      ArrayXd::LinSpaced(n_internal_samples_, bids.z_sets[0].minCoeff(),
                         bids.z_sets[0].maxCoeff());
  ArrayXd bidy_range =
      ArrayXd::LinSpaced(n_internal_samples_, bids.z_sets[1].minCoeff(),
                         bids.z_sets[1].maxCoeff());
  ArrayXXd interp_bidsx = Interpolate2D(bids.xs, bids.ys, bids.z_sets[0],
                                        value_pdfs_[id].xs, value_pdfs_[id].ys);
  ArrayXXd interp_bidsy = Interpolate2D(bids.xs, bids.ys, bids.z_sets[1],
                                        value_pdfs_[id].xs, value_pdfs_[id].ys);
  ArrayXXd cdf = TwoRandomVariableFunctionCDF(
      value_pdfs_[id].xs, value_pdfs_[id].ys, value_pdfs_[id].zs, interp_bidsx,
      interp_bidsy, bidx_range, bidy_range);
  bid_cdfs_[id] = {bidx_range, bidy_range, cdf};
  bid_margx_cdfs_[id] = {bidx_range, cdf.row(cdf.rows() - 1)};
  bid_margy_cdfs_[id] = {bidy_range, cdf.col(cdf.cols() - 1)};
}

float FirstPrice2D::GetFitness(const GridMulti& bids_in, int id) const {
  ArrayXd integrate_valsx =
      ArrayXd::LinSpaced((bids_in.xs.size() - 1) * 3, bids_in.xs(0),
                         bids_in.xs(bids_in.xs.size() - 1));
  ArrayXd integrate_valsy =
      ArrayXd::LinSpaced((bids_in.ys.size() - 1) * 3, bids_in.ys(0),
                         bids_in.ys(bids_in.ys.size() - 1));

  ArrayXXd valx_mesh = GetXMesh(integrate_valsx, integrate_valsy.size());
  ArrayXXd valy_mesh = GetYMesh(integrate_valsy, integrate_valsx.size());
  ArrayXXd bidsx = Interpolate2D(bids_in.xs, bids_in.ys, bids_in.z_sets[0],
                                 integrate_valsx, integrate_valsy);
  ArrayXXd bidsy = Interpolate2D(bids_in.xs, bids_in.ys, bids_in.z_sets[1],
                                 integrate_valsx, integrate_valsy);
  ArrayXXd win_both_probs = ArrayXXd::Ones(bidsx.rows(), bidsx.cols());
  ArrayXXd win_x_probs = ArrayXXd::Ones(bidsx.rows(), bidsx.cols());
  ArrayXXd win_y_probs = ArrayXXd::Ones(bidsx.rows(), bidsx.cols());

  for (int j = 0; j < n_players_; ++j) {
    if (j != id) {
      win_both_probs *= Interpolate2DGrid(bid_cdfs_[j].xs, bid_cdfs_[j].ys,
                                          bid_cdfs_[j].zs, bidsx, bidsy);
      for (int row = 0; row < win_x_probs.rows(); ++row) {
        win_x_probs.row(row) *= Interpolate(
            bid_margx_cdfs_[j].xs, bid_margx_cdfs_[j].ys, bidsx.row(row));
        win_y_probs.row(row) *= Interpolate(bid_margy_cdfs_[j], bidsy.row(row));
      }
    }
  }

  win_x_probs -= win_both_probs;
  win_y_probs -= win_both_probs;
  ArrayXXd zeros = ArrayXXd::Zero(bidsx.rows(), bidsx.cols());
  ArrayXXd win_none_probs = 1 - win_both_probs - win_x_probs - win_y_probs;
  ArrayXXd utils =
      (valx_mesh.binaryExpr(valy_mesh, utility_funcs_[id]) - bidsx - bidsy) *
          win_both_probs.unaryExpr(prob_weight_funcs_[id]) +
      (valx_mesh.binaryExpr(zeros, utility_funcs_[id]) - bidsx) *
          win_x_probs.unaryExpr(prob_weight_funcs_[id]) +
      (zeros.binaryExpr(valy_mesh, utility_funcs_[id]) - bidsy) *
          win_y_probs.unaryExpr(prob_weight_funcs_[id]) +
      utility_funcs_[id](0, 0) *
          (win_none_probs).unaryExpr(prob_weight_funcs_[id]);
  ArrayXXd likelihoods =
      Interpolate2D(value_pdfs_[id], integrate_valsx, integrate_valsy);
  return Areas2D(integrate_valsx, integrate_valsy, utils * likelihoods).sum();
}

float FirstPrice2D::GetExpectedRevenue(const GridMulti& bids_in, int id) const {
  ArrayXd integrate_valsx =
      ArrayXd::LinSpaced((bids_in.xs.size() - 1) * 3, bids_in.xs(0),
                         bids_in.xs(bids_in.xs.size() - 1));
  ArrayXd integrate_valsy =
      ArrayXd::LinSpaced((bids_in.ys.size() - 1) * 3, bids_in.ys(0),
                         bids_in.ys(bids_in.ys.size() - 1));

  ArrayXXd valx_mesh = GetXMesh(integrate_valsx, integrate_valsy.size());
  ArrayXXd valy_mesh = GetYMesh(integrate_valsy, integrate_valsx.size());
  ArrayXXd bidsx = Interpolate2D(bids_in.xs, bids_in.ys, bids_in.z_sets[0],
                                 integrate_valsx, integrate_valsy);
  ArrayXXd bidsy = Interpolate2D(bids_in.xs, bids_in.ys, bids_in.z_sets[1],
                                 integrate_valsx, integrate_valsy);
  ArrayXXd win_both_probs = ArrayXXd::Ones(bidsx.rows(), bidsx.cols());
  ArrayXXd win_x_probs = ArrayXXd::Ones(bidsx.rows(), bidsx.cols());
  ArrayXXd win_y_probs = ArrayXXd::Ones(bidsx.rows(), bidsx.cols());

  for (int j = 0; j < n_players_; ++j) {
    if (j != id) {
      win_both_probs *= Interpolate2DGrid(bid_cdfs_[j].xs, bid_cdfs_[j].ys,
                                          bid_cdfs_[j].zs, bidsx, bidsy);
      for (int row = 0; row < win_x_probs.rows(); ++row) {
        win_x_probs.row(row) *= Interpolate(
            bid_margx_cdfs_[j].xs, bid_margx_cdfs_[j].ys, bidsx.row(row));
        win_y_probs.row(row) *= Interpolate(bid_margy_cdfs_[j], bidsy.row(row));
      }
    }
  }

  win_x_probs -= win_both_probs;
  win_y_probs -= win_both_probs;
  ArrayXXd zeros = ArrayXXd::Zero(bidsx.rows(), bidsx.cols());
  ArrayXXd win_none_probs = 1 - win_both_probs - win_x_probs - win_y_probs;
  ArrayXXd utils =
      (bidsx + bidsy) * win_both_probs.unaryExpr(prob_weight_funcs_[id]) +
      (bidsx)*win_x_probs.unaryExpr(prob_weight_funcs_[id]) +
      (bidsy)*win_y_probs.unaryExpr(prob_weight_funcs_[id]) +
      0 * (win_none_probs).unaryExpr(prob_weight_funcs_[id]);
  ArrayXXd likelihoods =
      Interpolate2D(value_pdfs_[id], integrate_valsx, integrate_valsy);
  return Areas2D(integrate_valsx, integrate_valsy, utils * likelihoods).sum();
}

float FirstPrice2D::GetExpectedValue(const GridMulti& bids_in, int id) const {
  ArrayXd integrate_valsx =
      ArrayXd::LinSpaced((bids_in.xs.size() - 1) * 3, bids_in.xs(0),
                         bids_in.xs(bids_in.xs.size() - 1));
  ArrayXd integrate_valsy =
      ArrayXd::LinSpaced((bids_in.ys.size() - 1) * 3, bids_in.ys(0),
                         bids_in.ys(bids_in.ys.size() - 1));

  ArrayXXd valx_mesh = GetXMesh(integrate_valsx, integrate_valsy.size());
  ArrayXXd valy_mesh = GetYMesh(integrate_valsy, integrate_valsx.size());
  ArrayXXd bidsx = Interpolate2D(bids_in.xs, bids_in.ys, bids_in.z_sets[0],
                                 integrate_valsx, integrate_valsy);
  ArrayXXd bidsy = Interpolate2D(bids_in.xs, bids_in.ys, bids_in.z_sets[1],
                                 integrate_valsx, integrate_valsy);
  ArrayXXd win_both_probs = ArrayXXd::Ones(bidsx.rows(), bidsx.cols());
  ArrayXXd win_x_probs = ArrayXXd::Ones(bidsx.rows(), bidsx.cols());
  ArrayXXd win_y_probs = ArrayXXd::Ones(bidsx.rows(), bidsx.cols());

  for (int j = 0; j < n_players_; ++j) {
    if (j != id) {
      win_both_probs *= Interpolate2DGrid(bid_cdfs_[j].xs, bid_cdfs_[j].ys,
                                          bid_cdfs_[j].zs, bidsx, bidsy);
      for (int row = 0; row < win_x_probs.rows(); ++row) {
        win_x_probs.row(row) *= Interpolate(
            bid_margx_cdfs_[j].xs, bid_margx_cdfs_[j].ys, bidsx.row(row));
        win_y_probs.row(row) *= Interpolate(bid_margy_cdfs_[j], bidsy.row(row));
      }
    }
  }

  win_x_probs -= win_both_probs;
  win_y_probs -= win_both_probs;
  ArrayXXd zeros = ArrayXXd::Zero(bidsx.rows(), bidsx.cols());
  ArrayXXd win_none_probs = 1 - win_both_probs - win_x_probs - win_y_probs;
  ArrayXXd utils = (valx_mesh.binaryExpr(valy_mesh, utility_funcs_[id])) *
                       win_both_probs.unaryExpr(prob_weight_funcs_[id]) +
                   (valx_mesh.binaryExpr(zeros, utility_funcs_[id])) *
                       win_x_probs.unaryExpr(prob_weight_funcs_[id]) +
                   (zeros.binaryExpr(valy_mesh, utility_funcs_[id])) *
                       win_y_probs.unaryExpr(prob_weight_funcs_[id]) +
                   utility_funcs_[id](0, 0) *
                       (win_none_probs).unaryExpr(prob_weight_funcs_[id]);
  ArrayXXd likelihoods =
      Interpolate2D(value_pdfs_[id], integrate_valsx, integrate_valsy);
  return Areas2D(integrate_valsx, integrate_valsy, utils * likelihoods).sum();
}

}  // namespace auctions
