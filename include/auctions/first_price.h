#ifndef _AUCTIONS_FIRST_PRICE_H_
#define _AUCTIONS_FIRST_PRICE_H_

#include <functional>
#include <vector>

#include "numericaldists/distribution.h"

namespace auctions {

class FirstPrice {
 public:
  FirstPrice(std::vector<numericaldists::Distribution> value_dists);
  float GetValue(float value, float bid) const;
  void AcceptStrategy(std::function<float(float)> bid, int id);
  float GetFitness(const std::function<float(float)>& bid_func,
                    int id) const;

 private:
  float GetIntegrand(const std::function<float(float)>& bid_func, int id,
                      float value) const;
  std::vector<std::function<float(float)>> bid_funcs_;
  std::vector<std::function<float(float)>> bid_cdfs_;
  std::vector<numericaldists::Distribution> value_dists_;
  int n_players_;
};

}  // namespace auctions

#endif  // _AUCTIONS_FIRST_PRICE_H_
