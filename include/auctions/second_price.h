#ifndef AUCTIONS_SECOND_PRICE_H_
#define AUCTIONS_SECOND_PRICE_H_

#include <omp.h>
#include <functional>

#include "numericaldists/distribution.h"
#include "numericaldists/scatter.h"

#include <eigen3/Eigen/Core>

namespace auctions {

class SecondPrice {
 public:
  SecondPrice(std::vector<numericaldists::Distribution> value_dists,
              int n_internal_samples=10001);
  void AcceptStrategy(const numericaldists::Scatter& bids, int id);
  float GetFitness(const numericaldists::Scatter& bids, int id) const;
  float GetRevenue(const numericaldists::Scatter& bids, int id) const;
  float GetValue(const numericaldists::Scatter& bids, int id) const;
  std::vector<float> GetFitness(
      const std::vector<numericaldists::Scatter>& funcs, int id) const {
    if (!pre_calculated_) {
      Precalculate();
    }

    std::vector<float> fits(funcs.size(), 0.0);
#pragma omp parallel for
    for (int i = 0; i < funcs.size(); ++i) {
      fits[i] = GetFitness(funcs[i], id);
    }
    return fits;
  }

 private:
  void Precalculate() const;
  Eigen::ArrayXd internal_values_;
  Eigen::ArrayXd internal_bids_;
  std::vector<Eigen::ArrayXd> value_pdfs_;
  std::vector<Eigen::ArrayXd> bid_cdfs_;
  mutable std::vector<Eigen::ArrayXd> other_highest_cdfs_;
  mutable std::vector<Eigen::ArrayXd> exp_value_funcs_;
  std::vector<std::function<double(double)>> utility_funcs_;
  std::vector<std::function<double(double)>> prob_weight_funcs_;
  int n_players_;
  int n_internal_samples_;
  mutable bool pre_calculated_;
};

}  // namespace auctions

#endif  // AUCTIONS_SECOND_PRICE_H_
