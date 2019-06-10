#ifndef AUCTIONS_FIRST_PRICE_H_
#define AUCTIONS_FIRST_PRICE_H_

#include <functional>
#include <vector>

#include <iostream>

#include "numericaldists/distribution.h"
#include "numericaldists/function_ops.h"
#include "numericaldists/interval.h"
#include "numericaldists/scatter.h"

#include <eigen3/Eigen/Core>

namespace auctions {

class FirstPrice {
 public:
  FirstPrice(const std::vector<numericaldists::Distribution>& value_dists,
             int n_internal_samples = 10001);
  void AcceptStrategy(const numericaldists::Scatter& bids, int id);
  float GetFitness(const numericaldists::Scatter& bid_func, int id) const;
  float GetRevenue(const numericaldists::Scatter& bid_func, int id) const;
  float GetValue(const numericaldists::Scatter& bid_func, int id) const;

 private:
  float GetIntegrand(const Eigen::ArrayXd& bid_func, int id, float value) const;
  std::vector<numericaldists::Scatter> value_pdfs_;
  std::vector<numericaldists::Scatter> bid_cdfs_;
  std::vector<std::function<double(double)>> utility_funcs_;
  std::vector<std::function<double(double)>> prob_weight_funcs_;
  int n_players_;
  int n_internal_samples_;
  double max_bid_;
};

}  // namespace auctions

#endif  // AUCTIONS_FIRST_PRICE_H_
