#ifndef _NUMERICALDISTS_DISTRIBUTION_OPS_H_
#define _NUMERICALDISTS_DISTRIBUTION_OPS_H_

#include "numericaldists/bounds.h"
#include "numericaldists/distribution.h"
#include "numericaldists/interval.h"
#include "numericaldists/piecewise_linear.h"
#include "numericaldists/uneven_piecewise_linear.h"

namespace numericaldists {

UnevenPiecewiseLinear ApproximateRandomVariableFunctionCDF(
    const Distribution& dist, const std::function<float(float)>& func,
    int n_samples = 1001);

PiecewiseLinear ApproximateExpectedValueFunction(
    std::function<float(float)> pdf, std::function<float(float)> cdf,
    Interval interval, int n_samples = 1001);

PiecewiseLinear ApproximatePDFExpectedValueFunction(
    std::function<float(float)> pdf, Interval interval, int n_samples = 1001);

PiecewiseLinear ApproximateCDFExpectedValueFunction(
    std::function<float(float)> cdf, Interval interval, int n_samples = 1001);

PiecewiseLinear ApproximateRandomVariableFunctionCDF(
    const std::function<float(float, float)>& pdf,
    const std::function<float(float, float)>& func, Interval x_int,
    Interval y_int, Interval f_int, int n_samples = 1001);

PiecewiseLinear ApproximateRandomVariableFunctionCDF(
    const std::function<std::vector<std::vector<float>>(
        std::vector<float>, std::vector<float>)>& pdf,
    const std::function<std::vector<std::vector<float>>(
        std::vector<float>, std::vector<float>)>& func,
    Interval x_int, Interval y_int, Interval f_int, int n_samples = 1001);

std::function<float(float, float)> MultivariatePDFDomainTransform(
    const std::function<float(float, float)>& func,
    const std::function<float(float, float)>& x_trans,
    const std::function<float(float, float)>& y_trans,
    const std::function<float(float, float)>& jacobian, Interval x_int,
    Interval y_int);

// PiecewiseLinear ApproximateRandomVariableFunctionCDF(
//     const Distribution& dist, const std::function<float(float, float)>& func,
//     int n_samples = 1001);

}  // namespace numericaldists

#endif  // _NUMERICALDISTS_DISTRIBUTION_OPS_H_
