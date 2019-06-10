#ifndef NUMERICALDISTS_FUNCTION_OPS_H_
#define NUMERICALDISTS_FUNCTION_OPS_H_

#include "numericaldists/grid.h"
#include "numericaldists/interval.h"
#include "numericaldists/scatter.h"

#include <eigen3/Eigen/Dense>

namespace numericaldists {

Eigen::ArrayXd OverlapMeanPool(const Eigen::ArrayXd& ys);
Eigen::ArrayXXd OverlapMeanPool(const Eigen::ArrayXXd& zs);
Eigen::ArrayXd Differences(const Eigen::ArrayXd& xs);

Eigen::ArrayXd Interpolate(const Eigen::ArrayXd& xs, const Eigen::ArrayXd& ys,
                           const Eigen::ArrayXd& new_xs);

Eigen::ArrayXd Interpolate(const Scatter& points, const Eigen::ArrayXd& new_xs);
void Interpolate(const Scatter& points, Scatter& out);

std::vector<Eigen::ArrayXd> Interpolate(const Eigen::ArrayXd& xs,
                                        const std::vector<Eigen::ArrayXd>& ys,
                                        const Eigen::ArrayXd& new_xs);

Eigen::ArrayXd InterpolateUneven(const Eigen::ArrayXd& xs,
                                 const Eigen::ArrayXd& ys,
                                 const Eigen::ArrayXd& new_xs);
Eigen::ArrayXd InterpolateUneven(const Scatter& points,
                                 const Eigen::ArrayXd& new_xs);
void InterpolateUneven(const Scatter& points, Scatter& out);

std::vector<Eigen::ArrayXd> InterpolateUneven(
    const Eigen::ArrayXd& xs, const std::vector<Eigen::ArrayXd>& ys,
    const Eigen::ArrayXd& new_xs);

Eigen::ArrayXXd Interpolate2D(const Eigen::ArrayXd& xs,
                              const Eigen::ArrayXd& ys,
                              const Eigen::ArrayXXd& zs,
                              const Eigen::ArrayXd& new_xs,
                              const Eigen::ArrayXd& new_ys);
Eigen::ArrayXXd Interpolate2DGrid(const Eigen::ArrayXd& xs,
                              const Eigen::ArrayXd& ys,
                              const Eigen::ArrayXXd& zs,
                              const Eigen::ArrayXXd& new_xs,
                              const Eigen::ArrayXXd& new_ys);
Eigen::ArrayXXd Interpolate2D(const Grid& grid, const Eigen::ArrayXd& new_xs,
                              const Eigen::ArrayXd& new_ys);
void Interpolate2D(const Grid& grid, Grid& out);

Eigen::ArrayXXd Interpolate2DUneven(const Eigen::ArrayXd& xs,
                                    const Eigen::ArrayXd& ys,
                                    const Eigen::ArrayXXd& zs,
                                    const Eigen::ArrayXd& new_xs,
                                    const Eigen::ArrayXd& new_ys);
Eigen::ArrayXXd Interpolate2DUneven(const Grid& grid,
                                    const Eigen::ArrayXd& new_xs,
                                    const Eigen::ArrayXd& new_ys);
void Interpolate2DUneven(const Grid& grid, Grid& out);

std::vector<Eigen::ArrayXXd> Interpolate2D(
    const Eigen::ArrayXd& xs, const Eigen::ArrayXd& ys,
    const std::vector<Eigen::ArrayXXd>& zs, const Eigen::ArrayXd& new_xs,
    const Eigen::ArrayXd& new_ys);

std::vector<Eigen::ArrayXXd> Interpolate2DUneven(
    const Eigen::ArrayXd& xs, const Eigen::ArrayXd& ys,
    const std::vector<Eigen::ArrayXXd>& zs, const Eigen::ArrayXd& new_xs,
    const Eigen::ArrayXd& new_ys);

std::pair<Eigen::ArrayXd, Eigen::ArrayXd> InverseUneven(
    const Eigen::ArrayXd& xs, const Eigen::ArrayXd& ys);
std::pair<Eigen::ArrayXd, Eigen::ArrayXd> Inverse(const Eigen::ArrayXd& xs,
                                                  const Eigen::ArrayXd& ys,
                                                  Interval new_int);

Eigen::ArrayXd Derivative(const Eigen::ArrayXd& xs, const Eigen::ArrayXd& ys);
Eigen::ArrayXd Areas(const Eigen::ArrayXd& xs, const Eigen::ArrayXd& ys);
Eigen::ArrayXXd Areas2D(const Eigen::ArrayXd& xs, const Eigen::ArrayXd& ys,
                        const Eigen::ArrayXXd& zs);

Eigen::ArrayXd IntegralBelow(const Eigen::ArrayXd& xs,
                             const Eigen::ArrayXd& ys);
Eigen::ArrayXXd IntegralBelow2D(const Eigen::ArrayXd& xs,
                                const Eigen::ArrayXd& ys,
                                const Eigen::ArrayXXd& zs);
Eigen::ArrayXd IntegralAbove(const Eigen::ArrayXd& xs,
                             const Eigen::ArrayXd& ys);

}  // namespace numericaldists

#endif  // NUMERICALDISTS_BID_FUNCTION_OPS_H_
