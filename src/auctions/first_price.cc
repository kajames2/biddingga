#include "auctions/first_price.h"

#include <algorithm>

#include <boost/math/distributions/uniform.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include "numericaldists/distribution_ops.h"
#include "numericaldists/function_ops.h"
#include "numericaldists/order_statistic_ops.h"

using namespace boost::math;
using namespace numericaldists;

namespace auctions {

FirstPrice::FirstPrice(std::vector<Distribution> value_dists)
    : bid_funcs_(value_dists.size()),
      bid_cdfs_(value_dists.size()),
      value_dists_(value_dists),
      n_players_(value_dists.size()) {}

float FirstPrice::GetValue(float value, float bid) const { return value - bid; }

void FirstPrice::AcceptStrategy(std::function<float(float)> bid_func, int id) {
  bid_funcs_[id] = bid_func;
  bid_cdfs_[id] = ResampleFunction(
      ApproximateRandomVariableFunctionCDF(value_dists_[id], bid_func),
      Interval{0, upper(value_dists_)});
}

float FirstPrice::GetFitness(const std::function<float(float)>& bid_func,
                             int id) const {
  float exp_profit = quadrature::gauss_kronrod<float, 61>::integrate(
      [this, &bid_func, id](float value) {
        return GetIntegrand(bid_func, id, value);
      },
      lower(value_dists_[id]), upper(value_dists_[id]), 0, 0);
  return exp_profit;
}

float FirstPrice::GetIntegrand(const std::function<float(float)>& bid_func,
                               int id, float value) const {
  float bid = bid_func(value);
  float integrand = value - bid;
  integrand *= pdf(value_dists_[id], value);
  for (int j = 0; j < n_players_; ++j) {
    if (j != id) {
      integrand *= bid_cdfs_[j](bid);
    }
  }
  return integrand;
}

}  // namespace auctions
