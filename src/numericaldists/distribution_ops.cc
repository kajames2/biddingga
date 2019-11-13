#include <algorithm>
#include <functional>
#include <numeric>

#include "numericaldists/distribution.h"
#include "numericaldists/function_ops.h"

#include <cmath>
#include <eigen3/Eigen/Core>

using namespace Eigen;

namespace numericaldists {

ArrayXXd GetXMesh(const ArrayXd& xs, int y_size) {
  return xs.matrix().transpose().replicate(y_size, 1);
}

ArrayXXd GetYMesh(const ArrayXd& ys, int x_size) {
  return ys.matrix().replicate(1, x_size);
}

ArrayXd CDF(const ArrayXd& xs, const ArrayXd& pdf) {
  ArrayXd cdf = IntegralBelow(xs, pdf);
  cdf /= cdf(xs.size() - 1);
  return cdf;
}

ArrayXd PDF(const ArrayXd& xs, const ArrayXd& cdf) {
  return Derivative(xs, cdf);
}

ArrayXXd CDF2D(const ArrayXd& xs, const ArrayXd& ys, const ArrayXXd& pdf) {
  ArrayXXd cdf = IntegralBelow2D(xs, ys, pdf);
  cdf /= cdf(ys.size() - 1, xs.size() - 1);
  return cdf;
}

// This has pretty bad precision...
ArrayXXd PDF2D(const ArrayXd& xs, const ArrayXd& ys, const ArrayXXd& cdf) {
  ArrayXXd areas = ArrayXXd::Ones(cdf.rows(), cdf.cols());
  areas = Areas2D(xs, ys, areas);
  ArrayXXd pdf = cdf;
  for (int i = cdf.cols() - 1; i > 0; --i) {
    pdf.col(i) -= pdf.col(i - 1);
  }
  for (int i = cdf.rows() - 1; i > 0; --i) {
    pdf.row(i) -= pdf.row(i - 1);
  }
  pdf.bottomRightCorner(pdf.rows() - 1, pdf.cols() - 1) /=
      areas.bottomRightCorner(areas.rows() - 1, areas.cols() - 1);
  ArrayXXd centered_pdf = pdf.bottomRightCorner(pdf.rows() - 1, pdf.cols() - 1);
  ArrayXd centered_xs = OverlapMeanPool(xs);
  ArrayXd centered_ys = OverlapMeanPool(ys);
  pdf = Interpolate2DUneven(centered_xs, centered_ys, centered_pdf, xs, ys);
  pdf /= (pdf * areas).sum();
  return pdf;
}

// Esimates the CDF for a function of a random variable.
// new_xs must be linearly spaced.
ArrayXd RandomVariableFunctionCDF(const ArrayXd& xs, const ArrayXd& pdf,
                                  const ArrayXd& func_vals,
                                  const ArrayXd& new_xs) {
  const int size = xs.size();
  const int new_size = new_xs.size();
  const ArrayXd probs = Areas(xs, pdf).tail(size - 1);
  const ArrayXd f_vals = OverlapMeanPool(func_vals);
  const double spacing =
      (new_xs(new_xs.size() - 1) - new_xs(0)) / (new_size - 1);
  const ArrayXi inds = ((f_vals - new_xs(0)) / spacing)
                           .max(0)
                           .min(new_xs.size() - 1)
                           .ceil()
                           .cast<int>();
  ArrayXd f_pdf = ArrayXd::Zero(new_size);
  for (int i = 0; i < size - 1; ++i) {
    f_pdf[inds(i)] += probs(i);
  }
  std::partial_sum(f_pdf.data(), f_pdf.data() + new_size, f_pdf.data());
  f_pdf /= f_pdf(f_pdf.size() - 1);
  return f_pdf;
}

ArrayXd ExpectedValueFunction(const ArrayXd& xs, const ArrayXd& pdf,
                              const ArrayXd& cdf) {
  ArrayXd integrals = IntegralBelow(xs, pdf * xs);
  integrals /= cdf;
  integrals = integrals.min(xs).max(xs(0));
  return integrals;
}

ArrayXd ExpectedValueFunctionPDF(const ArrayXd& xs, const ArrayXd& pdf) {
  return ExpectedValueFunction(xs, pdf, CDF(xs, pdf));
}

ArrayXd ExpectedValueFunctionCDF(const ArrayXd& xs, const ArrayXd& cdf) {
  return ExpectedValueFunction(xs, PDF(xs, cdf), cdf);
}

ArrayXd ConditionalXPDF(const ArrayXd& xs, const ArrayXd& ys,
                        const ArrayXXd& pdf, float y) {
  ArrayXd y_arr(1);
  y_arr << y;
  ArrayXd pdf_slice = Interpolate2D(xs, ys, pdf, xs, y_arr);
  pdf_slice /= pdf_slice.sum();
  return pdf_slice;
}

ArrayXd ConditionalYPDF(const ArrayXd& xs, const ArrayXd& ys,
                        const ArrayXXd& pdf, float x) {
  ArrayXd x_arr(1);
  x_arr << x;
  ArrayXd pdf_slice = Interpolate2D(xs, ys, pdf, x_arr, ys);
  pdf_slice /= pdf_slice.sum();
  return pdf_slice;
}

ArrayXXd JointPDFIndependent(const ArrayXd& x_pdf, const ArrayXd& y_pdf) {
  return (y_pdf.matrix() * x_pdf.matrix().transpose()).array();
}

// Assume xs/ys/xs_f1/ys_f1 are linearly spaced
ArrayXXd TwoRandomVariableFunctionCDF(const ArrayXd& xs, const ArrayXd& ys,
                                      const ArrayXXd& joint_pdf,
                                      const ArrayXXd& f1, const ArrayXXd& f2,
                                      const ArrayXd& xs_f1,
                                      const ArrayXd& xs_f2) {
  assert(joint_pdf.size() == f1.size() && joint_pdf.size() == f2.size());

  ArrayXXd probs = Areas2D(xs, ys, joint_pdf)
                       .bottomRightCorner(ys.size() - 1, xs.size() - 1);
  ArrayXXd f1_vals = OverlapMeanPool(f1);
  double spacing_f1 = (xs_f1(xs_f1.size() - 1) - xs_f1(0)) / (xs_f1.size() - 1);
  ArrayXXi inds_f1 = ((f1_vals - xs_f1(0)) / spacing_f1)
                         .max(0)
                         .min(xs_f1.size() - 1)
                         .ceil()
                         .cast<int>();
  ArrayXXd f2_vals = OverlapMeanPool(f2);
  double spacing_f2 = (xs_f2(xs_f2.size() - 1) - xs_f2(0)) / (xs_f2.size() - 1);
  ArrayXXi inds_f2 = ((f2_vals - xs_f2(0)) / spacing_f2)
                         .max(0)
                         .min(xs_f2.size() - 1)
                         .ceil()
                         .cast<int>();
  ArrayXXd fs_cdf = ArrayXXd::Zero(xs_f2.size(), xs_f1.size());

  for (int row = 0; row < probs.rows(); ++row) {
    for (int col = 0; col < probs.cols(); ++col) {
      fs_cdf(inds_f2(row, col), inds_f1(row, col)) += probs(row, col);
    }
  }
  for (int row = 1; row < fs_cdf.rows(); ++row) {
    fs_cdf.row(row) += fs_cdf.row(row - 1);
  }
  for (int col = 1; col < fs_cdf.cols(); ++col) {
    fs_cdf.col(col) += fs_cdf.col(col - 1);
  }
  return fs_cdf;
}

ArrayXd RandomVariableFunctionCDF(const ArrayXd& xs, const ArrayXd& ys,
                                  const ArrayXXd& joint_pdf, const ArrayXXd& zs,
                                  const ArrayXd& new_xs) {
  int new_size = new_xs.size();
  ArrayXXd probs = Areas2D(xs, ys, joint_pdf)
                       .bottomRightCorner(ys.size() - 1, xs.size() - 1);
  ArrayXXd z_vals = OverlapMeanPool(zs);
  double z_spacing =
      (new_xs(new_xs.size() - 1) - new_xs(0)) / (new_xs.size() - 1);
  ArrayXXi inds = ((z_vals - new_xs(0)) / z_spacing)
                      .max(0)
                      .min(new_xs.size() - 1)
                      .ceil()
                      .cast<int>();

  ArrayXd f_pdf = ArrayXd::Zero(new_size);
  for (int row = 0; row < inds.rows() - 1; ++row) {
    for (int col = 0; col < inds.cols() - 1; ++col) {
      f_pdf[inds(row, col)] += probs(row, col);
    }
  }

  std::partial_sum(f_pdf.data(), f_pdf.data() + new_size, f_pdf.data());
  f_pdf /= f_pdf(f_pdf.size() - 1);
  return f_pdf;
}

// Assumes independence between pdf1 and pdf2
// ArrayXd RandomVariableBinaryOperationPDF(
//     const ArrayXd& xs, const ArrayXd& pdf1, const ArrayXd& pdf2,
//     const std::function<double(double, double)>& f, const ArrayXd& new_xs) {
//   double spacing = (new_xs(new_xs.size() - 1) - new_xs(0)) / new_xs.size();
//   ArrayXXd joint_pdf = JointPDFIndependent(pdf1, pdf2);
//   ArrayXXd x_mesh = GetXMesh(xs, xs.size());
//   ArrayXXd y_mesh = GetYMesh(xs, xs.size());
//   ArrayXXd f_mesh = x_mesh.binaryExpr(y_mesh, f);
//   //  f_mesh * joint_pdf
// }

ArrayXd ApplyConditional(const ArrayXd& xs, const ArrayXd& pdf,
                         const std::function<bool(double)>& cond) {
  ArrayXd pdf_cond =
      xs.unaryExpr([cond](double x) { return cond(x) ? 1.0 : 0.0; }) * pdf;
  pdf_cond /= Areas(xs, pdf_cond).sum();
  return pdf_cond;
}

ArrayXXd ApplyConditional(const ArrayXd& xs, const ArrayXd& ys,
                          const ArrayXXd& pdf,
                          const std::function<bool(double, double)>& cond) {
  MatrixXd x_mat = GetXMesh(xs, ys.size());
  MatrixXd y_mat = GetYMesh(ys, xs.size());
  ArrayXXd pdf_cond = x_mat
                          .binaryExpr(y_mat,
                                      [cond](double x, double y) {
                                        return cond(x, y) ? 0.0 : 1.0;
                                      })
                          .array() *
                      pdf;
  pdf_cond /= Areas2D(xs, ys, pdf_cond).sum();
  return pdf_cond;
}

// Requires xs equally spaced
ArrayXd MarginalX(const ArrayXd& xs, const ArrayXd& ys,
                  const ArrayXXd& joint_pdf) {
  ArrayXd x_likes = joint_pdf.colwise().sum();
  x_likes /= Areas(xs, x_likes).sum();
  return x_likes;
}

// Requires ys equally spaced
ArrayXd MarginalY(const ArrayXd& xs, const ArrayXd& ys,
                  const ArrayXXd& joint_pdf) {
  ArrayXd y_likes = joint_pdf.rowwise().sum();
  y_likes /= Areas(ys, y_likes).sum();
  return y_likes;
}

}  // namespace numericaldists
