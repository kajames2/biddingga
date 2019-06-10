#ifndef NUMERICALDISTS_DISTRIBUTION_OPS_H_
#define NUMERICALDISTS_DISTRIBUTION_OPS_H_

#include "numericaldists/interval.h"
#include "numericaldists/scatter.h"
#include "numericaldists/grid.h"

#include <eigen3/Eigen/Core>

namespace numericaldists {

Eigen::ArrayXXd GetXMesh(const Eigen::ArrayXd& xs, int y_size);
Eigen::ArrayXXd GetYMesh(const Eigen::ArrayXd& ys, int x_size);
Eigen::ArrayXd CDF(const Eigen::ArrayXd& xs, const Eigen::ArrayXd& pdf);
Eigen::ArrayXd PDF(const Eigen::ArrayXd& xs, const Eigen::ArrayXd& cdf);
Eigen::ArrayXXd CDF2D(const Eigen::ArrayXd& xs, const Eigen::ArrayXd& ys,
                      const Eigen::ArrayXXd& pdf);
Eigen::ArrayXXd PDF2D(const Eigen::ArrayXd& xs, const Eigen::ArrayXd& ys,
                      const Eigen::ArrayXXd& cdf);
Eigen::ArrayXXd JointPDFIndependent(const Eigen::ArrayXd& x_pdf,
                                    const Eigen::ArrayXd& y_pdf);

Eigen::ArrayXd RandomVariableFunctionCDF(const Eigen::ArrayXd& xs,
                                         const Eigen::ArrayXd& pdf,
                                         const Eigen::ArrayXd& func_vals,
                                         const Eigen::ArrayXd& new_xs);
Eigen::ArrayXd ExpectedValueFunction(const Eigen::ArrayXd& xs,
                                     const Eigen::ArrayXd& pdf,
                                     const Eigen::ArrayXd& cdf);
Eigen::ArrayXd ExpectedValueFunctionPDF(const Eigen::ArrayXd& xs,
                                        const Eigen::ArrayXd& pdf);
Eigen::ArrayXd ExpectedValueFunctionCDF(const Eigen::ArrayXd& xs,
                                        const Eigen::ArrayXd& cdf);
Eigen::ArrayXd ConditionalXPDF(const Eigen::ArrayXd& xs,
                               const Eigen::ArrayXd& ys,
                               const Eigen::ArrayXXd& pdf, float y);
Eigen::ArrayXd ConditionalYPDF(const Eigen::ArrayXd& xs,
                               const Eigen::ArrayXd& ys,
                               const Eigen::ArrayXXd& pdf, float x);
Eigen::ArrayXXd TwoRandomVariableFunctionCDF(const Eigen::ArrayXd& xs,
                                             const Eigen::ArrayXd& ys,
                                             const Eigen::ArrayXXd& joint_pdf,
                                             const Eigen::ArrayXXd& f1,
                                             const Eigen::ArrayXXd& f2,
                                             const Eigen::ArrayXd& xs_f1,
                                             const Eigen::ArrayXd& xs_f2);
Eigen::ArrayXd RandomVariableFunctionCDF(const Eigen::ArrayXd& xs,
                                         const Eigen::ArrayXd& ys,
                                         const Eigen::ArrayXXd& joint_pdf,
                                         const Eigen::ArrayXXd& zs,
                                         const Eigen::ArrayXd& new_xs);

Eigen::ArrayXd ApplyConditional(const Eigen::ArrayXd& xs,
                                const Eigen::ArrayXd& pdf,
                                const std::function<bool(double)>& cond);

Eigen::ArrayXXd ApplyConditional(const Eigen::ArrayXd& xs,
                                 const Eigen::ArrayXd& ys,
                                 const Eigen::ArrayXXd& pdf,
                                 const std::function<bool(double, double)>& cond);
Eigen::ArrayXd MarginalX(const Eigen::ArrayXd& xs, const Eigen::ArrayXd& ys,
                         const Eigen::ArrayXXd& joint_pdf);
Eigen::ArrayXd MarginalY(const Eigen::ArrayXd& xs, const Eigen::ArrayXd& ys,
                         const Eigen::ArrayXXd& joint_pdf);

}  // namespace numericaldists

#endif  // NUMERICALDISTS_DISTRIBUTION_OPS_H_
