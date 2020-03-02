#ifndef AUCTIONS_COMMON_VALUE_ENDPOINTS2_H_
#define AUCTIONS_COMMON_VALUE_ENDPOINTS2_H_

#include <omp.h>
#include <functional>
#include <vector>

#include "numericaldists/distribution.h"
#include "numericaldists/grid.h"
#include "numericaldists/interval.h"
#include "numericaldists/scatter.h"

#include <eigen3/Eigen/Core>

namespace auctions {

class CommonValueEndpoints2 {
 public:
  CommonValueEndpoints2(int n_bidders, numericaldists::Distribution value_dist,
                        numericaldists::Distribution error_dist,
                        int n_internal_samples = 2001,
                        int value_integration_samples = 501);
  void AcceptStrategy(numericaldists::Scatter bid_func, int id);
  void AcceptStrategy(float rel_bid_value, int id);
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
  Eigen::ArrayXd value_pdf_;
  std::vector<std::function<double(double)>> utility_funcs_;
  std::vector<std::function<double(double)>> prob_weight_funcs_;
  numericaldists::Distribution error_dist_;
};

}  // namespace auctions

#endif  // AUCTIONS_COMMON_VALUE_ENDPOINTS2_H_
