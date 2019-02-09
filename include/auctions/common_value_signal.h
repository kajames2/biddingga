#ifndef _AUCTIONS_COMMON_VALUE_SIGNAL_H_
#define _AUCTIONS_COMMON_VALUE_SIGNAL_H_

#include <functional>
#include <vector>
#include <omp.h>

#include "numericaldists/interval.h"
#include "numericaldists/piecewise_linear.h"
#include "numericaldists/bilerper.h"
#include "numericaldists/distribution.h"

namespace auctions {

class CommonValueSignal {
 public:
  CommonValueSignal(std::vector<int> n_draws, float epsilon,
                    numericaldists::Interval discount_int);
  void AcceptStrategy(numericaldists::PiecewiseLinear discount_func, int id);
  void AcceptStrategy(float discount_value, int id);
  float GetFitness(const std::function<float(float)>& discount_func,
                   int id) const;
  std::vector<float> GetFitness(
      const std::vector<numericaldists::PiecewiseLinear>& funcs, int id) const {
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
  float GetFitness(const float discount_value, int id) const;
  std::vector<float> GetFitness(
      const std::vector<float>& discounts, int id) const {
    if (!pre_calculated_) {
      Precalculate();
    }

    std::vector<float> fits(discounts.size(), 0.0);
    #pragma omp parallel for
    for (int i = 0; i < discounts.size(); ++i) {
      fits[i] = GetFitness(discounts[i], id);
    }
    return fits;
  }
  
 private:
  double GetIntegrand(const std::function<float(float)>& discount_func, int id,
                      float value) const;
  void Precalculate() const;

  int n_players_;
  mutable bool pre_calculated_;
  std::vector<int> n_draws_;
  std::vector<numericaldists::PiecewiseLinear> discount_funcs_;
  std::vector<float> one_draw_discounts_;
  mutable std::vector<std::function<float(float)>> others_bids_cdfs_;
  std::vector<numericaldists::Bilerper> value_dists_;
  std::vector<numericaldists::Distribution> one_draw_value_dists_;
  const numericaldists::Interval signal_int_;
  const numericaldists::Interval range_int_;
  const numericaldists::Interval bid_int_;
};

}  // namespace auctions

#endif  // _AUCTIONS_COMMON_VALUE_SIGNAL_H_
