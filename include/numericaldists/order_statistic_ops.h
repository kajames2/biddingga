#ifndef _NUMERICALDISTS_ORDER_STATISTIC_OPS_H_
#define _NUMERICALDISTS_ORDER_STATISTIC_OPS_H_

#include "numericaldists/bilerper.h"
#include "numericaldists/bounds.h"
#include "numericaldists/distribution.h"
#include "numericaldists/interval.h"
#include "numericaldists/piecewise_linear.h"

namespace numericaldists {

PiecewiseLinear ApproximateKthLowestOrderStatisticPDF(
    const std::function<float(float)>& cdf, Interval interval, int n_draws,
    int k, int n_samples = 1001);

PiecewiseLinear ApproximateKthLowestOrderStatisticCDF(
    const std::function<float(float)>& cdf, Interval interval, int n_draws,
    int k, int n_samples = 1001);

PiecewiseLinear ApproximateHighestOrderStatisticCDF(
    const std::function<float(float)> cdf, Interval interval, int n_draws,
    int n_samples = 1001);

PiecewiseLinear ApproximateKthLowestOrderStatisticPDF(
    const std::vector<std::function<float(float)>>& cdfs, Interval interval,
    int k, int n_samples = 2001);

PiecewiseLinear ApproximateKthLowestOrderStatisticCDF(
    const std::vector<std::function<float(float)>>& cdfs, Interval interval,
    int k, int n_samples = 2001);

PiecewiseLinear ApproximateHighestOrderStatisticCDF(
    const std::vector<std::function<float(float)>>& cdfs, Interval interval,
    int n_samples = 2001);

Bilerper ApproximateJointOrderStatisticPDF(const Distribution& dist,
                                           int n_draws, int j, int k,
                                           int n_samples = 101);

Bilerper ApproximateLowestHighestJointOrderStatisticPDF(
    const Distribution& dist, int n_draws, int n_samples = 101);

Bilerper ApproximateLowestHighestJointOrderStatisticCDF(
    const std::function<float(float)> cdf, Interval interval, int n_draws,
    int n_samples = 101);

}  // namespace numericaldists

#endif  // _NUMERICALDISTS_ORDER_STATISTIC_OPS_H_
