#ifndef AUCTIONS_FIRST_PRICE_REVERSE_H_
#define AUCTIONS_FIRST_PRICE_REVERSE_H_

#include <functional>
#include <vector>

#include <iostream>

#include "numericaldists/distribution.h"
#include "numericaldists/function_ops.h"
#include "numericaldists/interval.h"
#include "numericaldists/scatter.h"

#include <eigen3/Eigen/Core>

namespace auctions {

class FirstPriceReverse {
 public:
  FirstPriceReverse(const std::vector<numericaldists::Distribution>& cost_dists,
                    int n_internal_samples = 10001);
  void AcceptStrategy(const numericaldists::Scatter& bids, int id);
  float GetFitness(const numericaldists::Scatter& bid_func, int id) const;

 private:
  float GetIntegrand(const Eigen::ArrayXd& bid_func, int id, float value) const;
  std::vector<numericaldists::Scatter> cost_pdfs_;
  std::vector<numericaldists::Scatter> bid_cdfs_;
  std::vector<std::function<double(double)>> utility_funcs_;
  std::vector<std::function<double(double)>> prob_weight_funcs_;
  int n_players_;
  int n_internal_samples_;
};

}  // namespace auctions

#endif  // AUCTIONS_FIRST_PRICE_REVERSE_H_
