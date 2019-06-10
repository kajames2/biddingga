#include "numericaldists/function_ops.h"

#include <numeric>

#include "numericaldists/interval.h"

#include <eigen3/Eigen/Dense>

using namespace Eigen;

namespace numericaldists {

ArrayXd OverlapMeanPool(const ArrayXd& ys) {
  return (ys.head(ys.size() - 1) + ys.tail(ys.size() - 1)) / 2;
}

ArrayXXd OverlapMeanPool(const ArrayXXd& zs) {
  return (zs.topLeftCorner(zs.rows() - 1, zs.cols() - 1) +
          zs.topRightCorner(zs.rows() - 1, zs.cols() - 1) +
          zs.bottomLeftCorner(zs.rows() - 1, zs.cols() - 1) +
          zs.bottomRightCorner(zs.rows() - 1, zs.cols() - 1)) /
         4;
}

ArrayXd Differences(const ArrayXd& xs) {
  return xs.tail(xs.size() - 1) - xs.head(xs.size() - 1);
}

// Helper function for linear interpolation for an equally spaced grid.
// Calculates which points the interpolated values lie between in the grid and
// the 'alpha'---the fraction of the way to the next point.
std::pair<ArrayXi, ArrayXd> GetIndicesAlphas(Interval x_int,
                                             const ArrayXd& new_xs,
                                             int n_inds) {
  double span = GetSpan(x_int);
  ArrayXd inds_exact = ((new_xs - x_int.min) * (n_inds - 1) / span)
                           .max(0)
                           .min(n_inds - 1.000001);
  ArrayXi inds_lower = inds_exact.cast<int>();
  ArrayXd alphas = inds_exact - inds_lower.cast<double>();
  return {inds_lower, alphas};
}

std::pair<ArrayXi, ArrayXd> GetIndicesAlphasUneven(const ArrayXd& xs,
                                                   const ArrayXd& new_xs) {
  int start_i = 0;
  int end_i = new_xs.size() - 1;
  double x_min = xs(0);
  double x_max = xs(xs.size() - 1);
  ArrayXi inds(new_xs.size());
  ArrayXd alphas(new_xs.size());
  while (start_i <= end_i && new_xs(start_i) <= x_min) {
    inds(start_i) = 0;
    alphas(start_i) = 0;
    ++start_i;
  }
  while (end_i >= start_i && new_xs(end_i) >= x_max) {
    inds(end_i) = xs.size() - 2;
    alphas(end_i) = 1;
    --end_i;
  }
  ArrayXd diffs = Differences(xs);
  int n = 1;
  for (int i = start_i; i <= end_i; ++i) {
    while (new_xs(i) > xs(n)) {
      ++n;
    }
    alphas(i) = (new_xs(i) - xs(n - 1)) / diffs(n - 1);
    inds(i) = n - 1;
  }
  return {inds, alphas};
}

ArrayXd InterpolateWithAlphas(const ArrayXi& inds_l, const ArrayXd& alphas,
                              const ArrayXd& ys, const ArrayXd& new_xs) {
  ArrayXi inds_u = inds_l + 1;
  ArrayXd new_ys = (1 - alphas) * inds_l.unaryExpr([&ys](int i) {
    return ys(i);
  }) + alphas * inds_u.unaryExpr([&ys](int i) { return ys(i); });
  return new_ys;
}

ArrayXd Interpolate(Interval x_int, const ArrayXd& ys, const ArrayXd& new_xs) {
  auto [inds_l, alphas] = GetIndicesAlphas(x_int, new_xs, ys.size());
  return InterpolateWithAlphas(inds_l, alphas, ys, new_xs);
}

// When interpolating multiple functions using the same xs, only calculate the
// alphas once.
std::vector<ArrayXd> Interpolate(const ArrayXd& xs,
                                 const std::vector<ArrayXd>& ys,
                                 const ArrayXd& new_xs) {
  auto [inds_l, alphas] =
      GetIndicesAlphas(Interval{xs(0), xs(xs.size() - 1)}, new_xs, xs.size());
  std::vector<ArrayXd> new_ys;
  for (const auto& y : ys) {
    new_ys.push_back(InterpolateWithAlphas(inds_l, alphas, y, new_xs));
  }
  return new_ys;
}

ArrayXd Interpolate(const ArrayXd& xs, const ArrayXd& ys,
                    const ArrayXd& new_xs) {
  return Interpolate(Interval{xs(0), xs(xs.size() - 1)}, ys, new_xs);
}

ArrayXd InterpolateUneven(const ArrayXd& xs, const ArrayXd& ys,
                          const ArrayXd& new_xs) {
  auto [inds_l, alphas] = GetIndicesAlphasUneven(xs, new_xs);
  return InterpolateWithAlphas(inds_l, alphas, ys, new_xs);
}

std::vector<ArrayXd> InterpolateUneven(const ArrayXd& xs,
                                       const std::vector<ArrayXd>& ys,
                                       const ArrayXd& new_xs) {
  auto [inds_l, alphas] = GetIndicesAlphasUneven(xs, new_xs);
  std::vector<ArrayXd> new_ys;
  for (const auto& y : ys) {
    new_ys.push_back(InterpolateWithAlphas(inds_l, alphas, y, new_xs));
  }
  return new_ys;
}

ArrayXXd Interpolate2D(const ArrayXd& xs, const ArrayXd& ys, const ArrayXXd& zs,
                       const ArrayXd& new_xs, const ArrayXd& new_ys) {
  auto [inds_lx, alphasx] =
      GetIndicesAlphas({xs(0), xs(xs.size() - 1)}, new_xs, xs.size());
  auto [inds_ly, alphasy] =
      GetIndicesAlphas({ys(0), ys(ys.size() - 1)}, new_ys, ys.size());
  ArrayXi inds_uy = inds_ly + 1;
  ArrayXXd new_zs(new_ys.size(), new_xs.size());
  for (int i = 0; i < new_ys.size(); ++i) {
    ArrayXd interp_row =
        (1 - alphasy(i)) * zs.row(inds_ly(i)) + alphasy(i) * zs.row(inds_uy(i));
    new_zs.row(i) = InterpolateWithAlphas(inds_lx, alphasx, interp_row, new_xs);
  }
  return new_zs;
}

ArrayXXd Interpolate2DGrid(const ArrayXd& xs, const ArrayXd& ys,
                           const ArrayXXd& zs, const ArrayXXd& new_xs,
                           const ArrayXXd& new_ys) {
  double spanx = xs(xs.size() - 1) - xs(0);
  double spany = ys(ys.size() - 1) - ys(0);
  int nx_inds = ys.size();
  int ny_inds = xs.size();
  ArrayXXd indsx_exact =
      ((new_xs - xs(0)) * (nx_inds - 1) / spanx).max(0).min(nx_inds - 1.000001);
  ArrayXXd indsy_exact =
      ((new_ys - ys(0)) * (ny_inds - 1) / spany).max(0).min(ny_inds - 1.000001);
  ArrayXXi indsx_lower = indsx_exact.cast<int>();
  ArrayXXi indsy_lower = indsy_exact.cast<int>();
  ArrayXXd alphasx = indsx_exact - indsx_lower.cast<double>();
  ArrayXXd alphasy = indsy_exact - indsy_lower.cast<double>();

  ArrayXXi indsx_upper = indsx_lower + 1;
  ArrayXXi indsy_upper = indsy_lower + 1;
  ArrayXXd ll = indsy_lower.binaryExpr(indsx_lower,
                                       [zs](int x, int y) { return zs(y, x); });
  ArrayXXd lu = indsy_lower.binaryExpr(indsx_upper,
                                       [zs](int x, int y) { return zs(y, x); });
  ArrayXXd ul = indsy_upper.binaryExpr(indsx_lower,
                                       [zs](int x, int y) { return zs(y, x); });
  ArrayXXd uu = indsy_upper.binaryExpr(indsx_upper,
                                       [zs](int x, int y) { return zs(y, x); });
  ArrayXXd new_zs = (1 - alphasy) * (1 - alphasx) * ll +
                    (1 - alphasy) * alphasx * lu +
                    alphasy * (1 - alphasx) * ul + alphasy * alphasx * uu;
  return new_zs;
}

ArrayXXd Interpolate2DUneven(const ArrayXd& xs, const ArrayXd& ys,
                             const ArrayXXd& zs, const ArrayXd& new_xs,
                             const ArrayXd& new_ys) {
  auto [inds_lx, alphasx] = GetIndicesAlphasUneven(xs, new_xs);
  auto [inds_ly, alphasy] = GetIndicesAlphasUneven(ys, new_ys);
  ArrayXi inds_uy = inds_ly + 1;
  ArrayXXd new_zs(new_ys.size(), new_xs.size());
  for (int i = 0; i < new_ys.size(); ++i) {
    ArrayXd interp_row =
        (1 - alphasy(i)) * zs.row(inds_ly(i)) + alphasy(i) * zs.row(inds_uy(i));
    new_zs.row(i) = InterpolateWithAlphas(inds_lx, alphasx, interp_row, new_xs);
  }
  return new_zs;
}

std::vector<ArrayXXd> Interpolate2D(const ArrayXd& xs, const ArrayXd& ys,
                                    const std::vector<ArrayXXd>& zs,
                                    const ArrayXd& new_xs,
                                    const ArrayXd& new_ys) {
  auto [inds_lx, alphasx] =
      GetIndicesAlphas({xs(0), xs(xs.size() - 1)}, new_xs, xs.size());
  auto [inds_ly, alphasy] =
      GetIndicesAlphas({ys(0), ys(ys.size() - 1)}, new_ys, ys.size());
  ArrayXi inds_uy = inds_ly + 1;
  std::vector<ArrayXXd> new_zs;
  for (const auto& z : zs) {
    ArrayXXd new_z(new_ys.size(), new_xs.size());
    for (int i = 0; i < new_ys.size(); ++i) {
      ArrayXd interp_row =
          (1 - alphasy(i)) * z.row(inds_ly(i)) + alphasy(i) * z.row(inds_uy(i));
      new_z.row(i) =
          InterpolateWithAlphas(inds_lx, alphasx, interp_row, new_xs);
    }
    new_zs.push_back(new_z);
  }
  return new_zs;
}

std::vector<ArrayXXd> Interpolate2DUneven(const ArrayXd& xs, const ArrayXd& ys,
                                          const std::vector<ArrayXXd>& zs,
                                          const ArrayXd& new_xs,
                                          const ArrayXd& new_ys) {
  auto [inds_lx, alphasx] = GetIndicesAlphasUneven(xs, new_xs);
  auto [inds_ly, alphasy] = GetIndicesAlphasUneven(ys, new_ys);
  ArrayXi inds_uy = inds_ly + 1;
  std::vector<ArrayXXd> new_zs;
  for (const auto& z : zs) {
    ArrayXXd new_z(new_ys.size(), new_xs.size());
    for (int i = 0; i < new_ys.size(); ++i) {
      ArrayXd interp_row =
          (1 - alphasy(i)) * z.row(inds_ly(i)) + alphasy(i) * z.row(inds_uy(i));
      new_z.row(i) =
          InterpolateWithAlphas(inds_lx, alphasx, interp_row, new_xs);
    }
    new_zs.push_back(new_z);
  }
  return new_zs;
}

std::pair<ArrayXd, ArrayXd> InverseUneven(const ArrayXd& xs,
                                          const ArrayXd& ys) {
  ArrayXd new_xs = ys;
  ArrayXd new_ys = xs;
  if (new_xs(0) > new_xs(new_xs.size() - 1)) {
    new_xs.reverseInPlace();
    new_ys.reverseInPlace();
  }
  return {new_xs, new_ys};
}

std::pair<ArrayXd, ArrayXd> Inverse(const ArrayXd& xs, const ArrayXd& ys,
                                    Interval new_int) {
  auto [new_xs, new_ys] = InverseUneven(xs, ys);
  ArrayXd even_xs = ArrayXd::LinSpaced(ys.size(), new_int.min, new_int.max);
  new_ys = InterpolateUneven(new_xs, new_ys, even_xs);
  return {even_xs, new_ys};
}

ArrayXd Derivative(const ArrayXd& xs, const ArrayXd& ys) {
  int size = xs.size();
  ArrayXd x_mids = OverlapMeanPool(xs);
  ArrayXd derivs = Differences(ys) / Differences(xs);
  ArrayXd derivs_resample = Interpolate(x_mids, derivs, xs);
  return derivs_resample;
}

// Returns the area between each x with the previous x.  0 for first element.
// Trapazoidal method
ArrayXd Areas(const ArrayXd& xs, const ArrayXd& ys) {
  int size = xs.size();
  ArrayXd areas(size);
  areas(0) = 0;
  areas.tail(size - 1) = OverlapMeanPool(ys) * Differences(xs);
  return areas;
}

ArrayXXd Areas2D(const ArrayXd& xs, const ArrayXd& ys, const ArrayXXd& zs) {
  int x_size = xs.size();
  int y_size = ys.size();
  ArrayXXd areas = ArrayXXd::Zero(y_size, x_size);
  areas.block(1, 1, y_size - 1, x_size - 1) =
      OverlapMeanPool(zs) *
      (Differences(ys).matrix() * Differences(xs).matrix().transpose()).array();
  return areas;
}

ArrayXd IntegralBelow(const ArrayXd& xs, const ArrayXd& ys) {
  ArrayXd areas = Areas(xs, ys);
  std::partial_sum(areas.data(), areas.data() + xs.size(), areas.data());
  return areas;
}

ArrayXXd IntegralBelow2D(const ArrayXd& xs, const ArrayXd& ys,
                         const ArrayXXd& zs) {
  ArrayXXd areas = Areas2D(xs, ys, zs);
  for (int i = 1; i < ys.size(); ++i) {
    areas.row(i) += areas.row(i - 1);
  }
  for (int i = 1; i < xs.size(); ++i) {
    areas.col(i) += areas.col(i - 1);
  }
  return areas;
}

ArrayXd IntegralAbove(const ArrayXd& xs, const ArrayXd& ys) {
  ArrayXd int_below = IntegralBelow(xs, ys);
  ArrayXd int_above = int_below(xs.size() - 1) - int_below;
  return int_above;
}

Eigen::ArrayXd Interpolate(const Scatter& points,
                           const Eigen::ArrayXd& new_xs) {
  return Interpolate(points.xs, points.ys, new_xs);
}

void Interpolate(const Scatter& points, Scatter& out) {
  out.ys = Interpolate(points.xs, points.ys, out.xs);
}

Eigen::ArrayXd InterpolateUneven(const Scatter& points,
                                 const Eigen::ArrayXd& new_xs) {
  return Interpolate(points.xs, points.ys, new_xs);
}

void InterpolateUneven(const Scatter& points, Scatter& out) {
  out.ys = InterpolateUneven(points.xs, points.ys, out.xs);
}

Eigen::ArrayXXd Interpolate2D(const Grid& grid, const Eigen::ArrayXd& new_xs,
                              const Eigen::ArrayXd& new_ys) {
  return Interpolate2D(grid.xs, grid.ys, grid.zs, new_xs, new_ys);
}

void Interpolate2D(const Grid& grid, Grid& out) {
  out.zs = Interpolate2D(grid.xs, grid.ys, grid.zs, out.xs, out.ys);
}

}  // namespace numericaldists
