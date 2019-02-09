#ifndef _AUCTIONS_WINNER_PAY_H_
#define _AUCTIONS_WINNER_PAY_H_

#include <omp.h>
#include <functional>

#include "numericaldists/distribution.h"
#include "numericaldists/piecewise_linear.h"

namespace auctions {

class WinnerPay {
 public:
  WinnerPay(std::vector<numericaldists::Distribution> value_dists,
            int pay_nth_highest = 1);
  WinnerPay(std::vector<numericaldists::Distribution> value_dists,
            int pay_nth_highest, float max_bid);
  void AcceptStrategy(numericaldists::PiecewiseLinear bid, int id);
  float GetFitness(const numericaldists::PiecewiseLinear& bid_func,
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

 private:
  double GetIntegrand(const numericaldists::PiecewiseLinear& bid_func, int id,
                      float value) const;
  void Precalculate() const;
  std::vector<numericaldists::PiecewiseLinear> bid_funcs_;
  mutable std::vector<numericaldists::PiecewiseLinear>
      others_highest_order_stat_;
  mutable std::vector<numericaldists::PiecewiseLinear>
      others_nth_highest_exp_bid_;
  std::vector<numericaldists::Distribution> value_dists_;
  int n_players_;
  mutable bool pre_calculated_;
  int pay_nth_highest_;
  float max_bid_;
};

}  // namespace auctions

#endif  // _AUCTIONS_WINNER_PAY_H_
