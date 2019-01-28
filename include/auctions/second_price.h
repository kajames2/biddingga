#ifndef _AUCTIONS_SECOND_PRICE_H_
#define _AUCTIONS_SECOND_PRICE_H_

#include <functional>

#include "numericaldists/distribution.h"

namespace auctions {

class SecondPrice {
 public:
  SecondPrice(std::vector<numericaldists::Distribution> value_dists);
  void AcceptStrategy(std::function<float(float)> bid, int id);
  float GetFitness(const std::function<float(float)>& bid_func,
                    int id) const;

 private:
  double GetIntegrand(const std::function<float(float)>& bid_func, int id,
                      float value) const;
  void Precalculate() const;
  std::vector<std::function<float(float)>> bid_funcs_;
  mutable std::vector<std::function<float(float)>> order_stat_funcs_;
  mutable std::vector<std::function<float(float)>> exp_value_funcs_;
  std::vector<numericaldists::Distribution> value_dists_;
  int n_players_;
  mutable bool pre_calculated_;
};

}  // namespace auctions

#endif  // _AUCTIONS_SECOND_PRICE_H_
