#ifndef AUCTIONS_COMMON_VALUE_SIGNAL_H_
#define AUCTIONS_COMMON_VALUE_SIGNAL_H_

#include <omp.h>
#include <functional>
#include <vector>

#include "numericaldists/distribution.h"
#include "numericaldists/grid.h"
#include "numericaldists/interval.h"
#include "numericaldists/scatter.h"

#include <eigen3/Eigen/Core>

namespace auctions {

class CommonValueSignal {
 public:
  CommonValueSignal(std::vector<int> n_draws, float epsilon,
                    numericaldists::Interval rel_bid_int,
                    int n_internal_samples = 1001,
                    int mstar_integration_samples = 121);
  void AcceptStrategy(numericaldists::Scatter rel_bid_func, int id);
  void AcceptStrategy(float rel_bid_value, int id);
  float GetFitness(const numericaldists::Scatter& rel_bid_func, int id) const;
  std::vector<float> GetFitness(
      const std::vector<numericaldists::Scatter>& rel_bid_funcs,
      int id) const {
    if (!pre_calculated_) {
      Precalculate();
    }

    std::vector<float> fits(rel_bid_funcs.size(), 0.0);
    #pragma omp parallel for
    for (int i = 0; i < rel_bid_funcs.size(); ++i) {
      fits[i] = GetFitness(rel_bid_funcs[i], id);
    }
    return fits;
  }

  float GetFitness(const float rel_bid_value, int id) const;
  std::vector<float> GetFitness(const std::vector<float>& rel_bids,
                                int id) const {
    if (!pre_calculated_) {
      Precalculate();
    }

    std::vector<float> fits(rel_bids.size(), 0.0);
#pragma omp parallel for
    for (int i = 0; i < rel_bids.size(); ++i) {
      fits[i] = GetFitness(rel_bids[i], id);
    }
    return fits;
  }

 private:
  double GetIntegrand(const numericaldists::Scatter& rel_bid_func, int id,
                      float value) const;
  void Precalculate() const;

  int n_players_;
  mutable bool pre_calculated_;
  std::vector<int> n_draws_;
  std::vector<numericaldists::Scatter> rel_bid_funcs_;
  std::vector<float> one_draw_rel_bids_;
  mutable std::vector<Eigen::ArrayXd> others_bids_cdfs_;
  Eigen::ArrayXd internal_mstars_;
  Eigen::ArrayXd internal_precs_;
  Eigen::ArrayXd internal_bids_;
  Eigen::ArrayXd internal_signals_;
  std::vector<Eigen::ArrayXXd> value_dists_;
  Eigen::ArrayXd one_draw_pdf_;
  int n_internal_samples_;
  int mstar_integration_samples_;
  std::vector<std::function<double(double)>> utility_funcs_;
  std::vector<std::function<double(double)>> prob_weight_funcs_;
};

}  // namespace auctions

#endif  // AUCTIONS_COMMON_VALUE_SIGNAL_H_
