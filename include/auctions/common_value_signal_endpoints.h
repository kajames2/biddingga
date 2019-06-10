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
                             std::vector<int> n_draws,
                             float epsilon,
                             int n_internal_samples = 201,
                             int value_integration_samples = 201);
  void AcceptStrategy(numericaldists::Grid bid_func, int id);
  void AcceptStrategy(float rel_bid_value, int id);
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

 private:
  double GetIntegrand(const numericaldists::Scatter& rel_bid_func, int id,
                      float value) const;
  void Precalculate() const;

  int n_players_;
  mutable bool pre_calculated_;
  mutable std::vector<Eigen::ArrayXXd> others_bids_cdfs_;
  std::vector<numericaldists::Scatter> bid_funcs_;
  Eigen::ArrayXd internal_values_;
  Eigen::ArrayXd internal_signals_;
  Eigen::ArrayXd internal_bids_;
  Eigen::ArrayXXd value_pdf_;
  std::vector<std::function<double(double)>> utility_funcs_;
  std::vector<std::function<double(double)>> prob_weight_funcs_;
};

}  // namespace auctions

#endif  // AUCTIONS_COMMON_VALUE_SIGNAL_ENDPOINTS_H_
