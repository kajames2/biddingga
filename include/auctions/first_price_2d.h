#ifndef AUCTIONS_FIRST_PRICE_2D_H_
#define AUCTIONS_FIRST_PRICE_2D_H_

#include <functional>
#include <vector>

#include <iostream>

#include "numericaldists/distribution.h"
#include "numericaldists/function_ops.h"
#include "numericaldists/grid_multi.h"
#include "numericaldists/interval.h"

#include <eigen3/Eigen/Core>

namespace auctions {

class FirstPrice2D {
 public:
  FirstPrice2D(const std::vector<numericaldists::Distribution>& valuex_dists,
               const std::vector<numericaldists::Distribution>& valuey_dists,
               std::vector<std::function<double(double, double)>> utils,
               int n_internal_samples = 501);
  FirstPrice2D(const std::vector<numericaldists::Distribution>& valuex_dists,
               const std::vector<numericaldists::Distribution>& valuey_dists,
               int n_internal_samples = 501);

  void AcceptStrategy(const numericaldists::GridMulti& bids, int id);
  float GetFitness(const numericaldists::GridMulti& bid_func, int id) const;
  float GetExpectedRevenue(const numericaldists::GridMulti& bid_func,
                           int id) const;
  float GetExpectedValue(const numericaldists::GridMulti& bid_func,
                         int id) const;
 private:
  void Initialize(
      const std::vector<numericaldists::Distribution>& valuex_dists,
      const std::vector<numericaldists::Distribution>& valuey_dists);
  float GetIntegrand(const Eigen::ArrayXd& bid_func, int id, float value) const;
  std::vector<numericaldists::Grid> value_pdfs_;
  std::vector<numericaldists::Grid> bid_cdfs_;
  std::vector<std::function<double(double, double)>> utility_funcs_;
  std::vector<std::function<double(double)>> prob_weight_funcs_;
  std::vector<numericaldists::Scatter> bid_margx_cdfs_;
  std::vector<numericaldists::Scatter> bid_margy_cdfs_;
  int n_players_;
  int n_internal_samples_;
};

}  // namespace auctions

#endif  // AUCTIONS_FIRST_PRICE_2D_H_
