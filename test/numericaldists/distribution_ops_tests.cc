#include <gtest/gtest.h>

#include <vector>

#include <boost/math/distributions/normal.hpp>
#include "numericaldists/distribution.h"
#include "numericaldists/distribution_ops.h"
#include "numericaldists/function_ops.h"

#include <eigen3/Eigen/Core>

namespace gatests {

using namespace numericaldists;
using namespace Eigen;

class DistributionOpsTest : public ::testing::Test {
 public:
  DistributionOpsTest() {}
  ArrayXd xs = ArrayXd::LinSpaced(101, -1, 1);
  
  std::function<double(double)> sqr = [](double x) -> double {
    if (x < -1 || x > 1) { return 0; }
    return 1.5 * x * x;
  };
  
  std::function<double(double)> cub = [](double x) -> double {
    if (x < -1) { return 0; }
    else if (x > 1) { return 1;}
    else { return 0.5 * (1 + x * x * x); }
  };

  ArrayXd sqr_s = xs.unaryExpr(sqr);
  ArrayXd cub_s = xs.unaryExpr(cub);

 protected:
  virtual void SetUp() {}
};

TEST_F(DistributionOpsTest, CDF2DTest) {
  ArrayXd x = ArrayXd::LinSpaced(101, 0, 1);
  ArrayXd y = ArrayXd::LinSpaced(101, 0, 1);
  VectorXd xsqr = x.square().matrix();
  VectorXd ysqr = y.square().matrix();
  ArrayXXd zs = 3/2. * (xsqr.transpose().replicate(ysqr.size(), 1) +
                ysqr.replicate(1, xsqr.size()));
  ArrayXXd areas = CDF2D(x, y, zs);
  float epsilon = 0.001;
  EXPECT_NEAR(0, areas(0, 0), epsilon);
  EXPECT_NEAR(1, areas(100, 100), epsilon);
  EXPECT_NEAR(0, areas(100, 0), epsilon);
  EXPECT_NEAR(0, areas(0, 100), epsilon);
  EXPECT_NEAR(0.0625, areas(50, 50), epsilon);
  EXPECT_NEAR(0.3125, areas(100, 50), epsilon);
}

TEST_F(DistributionOpsTest, PDF2DTest) {
  ArrayXd x = ArrayXd::LinSpaced(101, 0, 1);
  ArrayXd y = ArrayXd::LinSpaced(101, 0, 1);
  VectorXd xsqr = x.square().matrix();
  VectorXd ysqr = y.square().matrix();
  ArrayXXd zs = 3/2. * (xsqr.transpose().replicate(ysqr.size(), 1) +
                ysqr.replicate(1, xsqr.size()));
  ArrayXXd pdf = PDF2D(x, y, CDF2D(x, y, zs));
  float epsilon = 0.05;
  EXPECT_NEAR(0, pdf(0, 0), epsilon);
  EXPECT_NEAR(2.94, pdf(99, 99), epsilon);
  EXPECT_NEAR(1.5, pdf(100, 0), epsilon);
  EXPECT_NEAR(1.5, pdf(0, 100), epsilon);
  EXPECT_NEAR(0.75, pdf(50, 50), epsilon);
  EXPECT_NEAR(1.875, pdf(100, 50), epsilon);
}

TEST_F(DistributionOpsTest, ExpectedValueFunctionPDFTest) {
  ArrayXd exp_vals = ExpectedValueFunctionPDF(xs, sqr_s);
  ArrayXd test_xs(3);
  test_xs << -1, 0, 1;
  ArrayXd test_ys = Interpolate(xs, exp_vals, test_xs);
  float epsilon = 0.001;
  EXPECT_NEAR(-1, test_ys(0), epsilon);
  EXPECT_NEAR(-3 / 4., test_ys(1), epsilon);
  EXPECT_NEAR(0, test_ys(2), epsilon);
}

TEST_F(DistributionOpsTest, CDFExpectedValueFunctionTest) {
  ArrayXd exp_vals = ExpectedValueFunctionCDF(xs, cub_s);
  ArrayXd test_xs(3);
  test_xs << -1, 0, 1;
  ArrayXd test_ys = Interpolate(xs, exp_vals, test_xs);
  float epsilon = 0.001;
  EXPECT_NEAR(-1, test_ys(0), epsilon);
  EXPECT_NEAR(-3 / 4., test_ys(1), epsilon);
  EXPECT_NEAR(0, test_ys(2), epsilon);
}

TEST_F(DistributionOpsTest, RandomVariableFunctionCDFTest) {
  Distribution dist = boost::math::normal_distribution<>(0, 1);
  ArrayXd xs = ArrayXd::LinSpaced(101, lower(dist), upper(dist));
  ArrayXd pdf_s = xs.unaryExpr([&dist](double x) -> double {return pdf(dist, x);});
  auto lin = [](double x) { return 2 * x; };
  ArrayXd lin_s = xs.unaryExpr(lin);
  ArrayXd new_xs = ArrayXd::LinSpaced(101, lin(lower(dist)), lin(upper(dist)));
  auto rv_cdf = RandomVariableFunctionCDF(xs, pdf_s, lin_s, new_xs);
  ArrayXd test_xs(3);
  test_xs << -4, 0, 2;
  ArrayXd test_ys = Interpolate(new_xs, rv_cdf, test_xs);
  float epsilon = 0.001;
  EXPECT_NEAR(0.022750, test_ys(0), epsilon);
  EXPECT_NEAR(1 / 2., test_ys(1), epsilon);
  EXPECT_NEAR(0.841345, test_ys(2), epsilon);
}

// TEST_F(DistributionOpsTest, RandomVariableFunctionCDF2DTest) {
//   auto pdf = [](float x, float y) { return x / 2 + 2 * y; };
//   auto func = [](float x, float y) { return 1 / 4. * x * x - 4 * y * y; };
//   auto out = RandomVariableFunctionCDF(pdf, func, {0, 2}, {0, 0.5}, {-1, 1});
//   float epsilon = 0.001;
//   EXPECT_NEAR(0, out(-2), epsilon);
//   EXPECT_NEAR(0, out(-1), epsilon);
//   EXPECT_NEAR(0.152369, out(-.5), epsilon);
//   EXPECT_NEAR(1 / 2., out(0), epsilon);
//   EXPECT_NEAR(1, out(1), epsilon);
//   EXPECT_NEAR(1, out(2), epsilon);
// }

// std::function<float(float, float)> MultivariatePDFDomainTransform(
//     const std::function<float(float, float)>& func,
//     const std::function<float(float, float)>& x_trans,
//     const std::function<float(float, float)>& y_trans,
//     const std::function<float(float, float)>& jacobian);

}  // namespace gatests
