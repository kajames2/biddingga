#ifndef _AUCTIONS_COMMON_VALUE_SIGNAL_H_
#define _AUCTIONS_COMMON_VALUE_SIGNAL_H_

#include <functional>
#include <vector>

#include "numericaldists/interval.h"

namespace auctions {

class CommonValueSignal {
 public:
  CommonValueSignal(std::vector<int> n_draws, float epsilon,
                    numericaldists::Interval discount_int);
  void AcceptStrategy(std::function<float(float)> discount_func, int id);
  float GetFitness(const std::function<float(float)>& discount_func,
                   int id) const;

 private:
  double GetIntegrand(const std::function<float(float)>& discount_func, int id,
                      float value) const;
  void Precalculate() const;

  int n_players_;
  mutable bool pre_calculated_;
  std::vector<int> n_draws_;
  float epsilon_;
  std::vector<std::function<float(float)>> discount_funcs_;
  mutable std::vector<std::function<float(float)>> others_bids_cdfs_;
  std::vector<std::function<float(float, float)>> value_dists_;
  const numericaldists::Interval signal_int_;
  const numericaldists::Interval range_int_;
  const numericaldists::Interval bid_int_;
};

}  // namespace auctions

#endif  // _AUCTIONS_COMMON_VALUE_SIGNAL_H_
