#ifndef AUCTIONS_COMMON_VALUE_SIGNAL_ENDPOINTS_H_
#define AUCTIONS_COMMON_VALUE_SIGNAL_ENDPOINTS_H_

#include <omp.h>
#include <functional>
#include <vector>

#include "numericaldists/distribution.h"
#include "numericaldists/grid.h"
#include "numericaldists/interval.h"
#include "numericaldists/scatter.h"

#include <eigen3/Eigen/Core>

namespace auctions {

class CommonValueSignalEndpoints {
 public:
  CommonValueSignalEndpoints(numericaldists::Distribution value_dist,
                             numericaldists::Distribution error_dist,
                             std::vector<int> n_draws,
                             int n_internal_samples = 201,
                             int value_integration_samples = 201);
  void AcceptStrategy(numericaldists::Grid bid_func, int id);
  void AcceptStrategy(numericaldists::Scatter bid, int id);
  float GetFitness(const numericaldists::Grid& bid_func, int id) const;
  std::vector<float> GetFitness(
      const std::vector<numericaldists::Grid>& bid_funcs, int id) const {
    if (!pre_calculated_) {
      Precalculate();
    }

    std::vector<float> fits(bid_funcs.size(), 0.0);
#pragma omp parallel for
    for (int i = 0; i < bid_funcs.size(); ++i) {
      fits[i] = GetFitness(bid_funcs[i], id);
    }
    return fits;
  }

  float GetFitness(const numericaldists::Scatter& bid_func, int id) const;
  std::vector<float> GetFitness(
      const std::vector<numericaldists::Scatter>& bid_funcs, int id) const {
    if (!pre_calculated_) {
      Precalculate();
    }

    std::vector<float> fits(bid_funcs.size(), 0.0);
#pragma omp parallel for
    for (int i = 0; i < bid_funcs.size(); ++i) {
      fits[i] = GetFitness(bid_funcs[i], id);
    }
    return fits;
  }

  Eigen::ArrayXXd GetSignalPDF(int id, int value_index,
                               Eigen::ArrayXd midpoints,
                               Eigen::ArrayXd precisions) const;

 private:
  double GetIntegrand(const numericaldists::Grid& bid_func, int id,
                      float value) const;
  void Precalculate() const;

  int n_players_;
  int n_internal_samples_;
  std::vector<int> n_draws_;
  mutable bool pre_calculated_;
  std::vector<numericaldists::Grid> bid_funcs_;
  std::vector<numericaldists::Scatter> one_draw_bids_;
  mutable std::vector<Eigen::ArrayXXd> others_bids_cdfs_;
  Eigen::ArrayXd internal_values_;
  Eigen::ArrayXd internal_midpoints_;
  Eigen::ArrayXd internal_precisions_;
  Eigen::ArrayXd internal_bids_;
  Eigen::ArrayXd value_pdf_;
  std::vector<numericaldists::Grid> relative_pdfs_;
  std::vector<std::function<double(double)>> utility_funcs_;
  std::vector<std::function<double(double)>> prob_weight_funcs_;
  numericaldists::Distribution error_dist_;
};

}  // namespace auctions

#endif  // AUCTIONS_COMMON_VALUE_SIGNAL_ENDPOINTS_H_
