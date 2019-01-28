#ifndef _NUMERICALDISTS_FUNCTION_OPS_H_
#define _NUMERICALDISTS_FUNCTION_OPS_H_

#include "numericaldists/bilerper.h"
#include "numericaldists/interval.h"
#include "numericaldists/piecewise_linear.h"
#include "numericaldists/uneven_piecewise_linear.h"

namespace numericaldists {

PiecewiseLinear ResampleFunction(const std::function<float(float)>& func,
                                 Interval interval, int n_samples = 10001);
PiecewiseLinear ResampleFunction(const UnevenPiecewiseLinear& func,
                                 Interval interval, int n_samples = 10001);

UnevenPiecewiseLinear ApproximateInverse(
    const std::function<float(float)>& func, Interval interval,
    int n_samples = 1001);

PiecewiseLinear ApproximateDerivative(const std::function<float(float)>& func,
                                      Interval interval, int n_samples = 1001);

std::vector<float> ApproximateAreas(std::function<float(float)> func,
                                    const std::vector<float>& x_samples,
                                    std::function<float(float)> integrand_func =
                                        [](float x) { return 1; });

PiecewiseLinear ApproximateIntegralBelow(
    std::function<float(float)> func, Interval interval, int n_samples = 1001,
    std::function<float(float)> integrand_func = [](float x) { return 1; });

PiecewiseLinear ApproximateIntegralAbove(
    std::function<float(float)> func, Interval interval, int n_samples = 1001,
    std::function<float(float)> integrand_func = [](float x) { return 1; });

Bilerper ResampleFunction2D(const std::function<float(float, float)>& func,
                            Interval x_int, Interval y_int,
                            int n_samples = 1001);

}  // namespace numericaldists

#endif  // _NUMERICALDISTS_BID_FUNCTION_OPS_H_
