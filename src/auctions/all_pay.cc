#include "auctions/all_pay.h"

#include <boost/math/distributions/uniform.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include "numericaldists/distribution_ops.h"
#include "numericaldists/function_ops.h"
#include "numericaldists/interval.h"
#include "numericaldists/order_statistic_ops.h"

using namespace boost::math;
using namespace numericaldists;

namespace auctions {

AllPay::AllPay(std::vector<float> values)
    : bid_cdf_funcs_(values.size()),
      bid_pdf_funcs_(values.size()),
      values_(values),
      n_players_(values.size()) {}

float AllPay::GetValue(float bid, float value) const { return value; }

void AllPay::AcceptStrategy(std::function<float(float)> cdf, int id) {
  bid_cdf_funcs_[id] = cdf;
  bid_pdf_funcs_[id] = ApproximateDerivative(cdf, Interval{0, values_[id]});
}

float ExpectedProfitTiesAtZero(
    const std::vector<std::function<float(float)>>& bid_cdfs, float value,
    int id) {
  int n_players = bid_cdfs.size();
  float prob_all_zero = 1;
  for (int j = 0; j < n_players; ++j) {
    prob_all_zero *= bid_cdfs[j](0);
  }
  return 1. / n_players * value * prob_all_zero;
}

float AllPay::GetFitness(const std::function<float(float)>& cdf_func,
                         int id) const {
  float exp_profit = quadrature::gauss_kronrod<float, 61>::integrate(
      [this, &cdf_func, id](float bid) {
        return GetIntegrand(cdf_func, id, bid);
      },
      0, values_[id], 0, 0);

  exp_profit += ExpectedProfitTiesAtZero(bid_cdf_funcs_, values_[id], id);
  float prob_no_other_val = 1;
  float prob_no_above = 1;
  for (int j = 0; j < n_players_; ++j) {
    if (j != id) {
      if (values_[j] == values_[id]) {
        prob_no_other_val *= bid_cdf_funcs_[j](values_[id]);
      } else {
        prob_no_above *= bid_cdf_funcs_[j](values_[id]);
      }
    }
  }
  float prob_val = (1 - cdf_func(values_[id]));
  exp_profit += (0 - values_[id]) * (1 - prob_no_above) * prob_val;
  exp_profit += ((1. / n_players_ * values_[id]) - values_[id]) *
                (1 - prob_no_other_val) * prob_no_above * prob_val;
  // exp_profit += (values_[id] - values_[id]) * prob_no_other_val *
  // prob_no_above * prob_val;
  return exp_profit;
}

float AllPay::GetIntegrand(const std::function<float(float)>& cdf_func, int id,
                           float bid) const {
  float density = bid_pdf_funcs_[id](bid);
  float prob_win = 1;
  for (int j = 0; j < n_players_; ++j) {
    if (j != id) {
      prob_win *= bid_cdf_funcs_[j](bid);
    }
  }
  return (prob_win * values_[id] - bid) * density;
}

}  // namespace auctions
