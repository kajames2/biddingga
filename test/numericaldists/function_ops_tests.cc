#include <gtest/gtest.h>

#include "numericaldists/function_ops.h"

#include <eigen3/Eigen/Core>

namespace gatests {

using namespace numericaldists;
using namespace Eigen;

ArrayXd xs_lofi = ArrayXd::LinSpaced(11, 0, 10);
ArrayXd sqr_lofi = xs_lofi * xs_lofi;
ArrayXd xs = ArrayXd::LinSpaced(10001, 0, 10);
ArrayXd horz = ArrayXd::Constant(10001, 5);
ArrayXd lin = 2 * xs;
ArrayXd sqr = xs * xs;
ArrayXd cube = xs * xs * xs;

class FunctionOpsTest : public ::testing::Test {
 public:
  FunctionOpsTest() {}

 protected:
  virtual void SetUp() {}
};

// Eigen::ArrayXd InterpolateUneven(const Eigen::ArrayXd& xs,
//                                  const Eigen::ArrayXd& ys,
//                                  const Eigen::ArrayXd& new_xs);

TEST_F(FunctionOpsTest, InterpolateTest) {
  ArrayXd ys = Interpolate(xs_lofi, sqr_lofi, ArrayXd::LinSpaced(101, 0, 10));
  float epsilon = 0.005;
  EXPECT_NEAR(0.1, ys(1), epsilon);
  EXPECT_NEAR(1, ys(10), epsilon);
  EXPECT_NEAR(2.2, ys(14), epsilon);
  EXPECT_NEAR(72.5, ys(85), epsilon);
}


TEST_F(FunctionOpsTest, InterpolateUnevenTest) {
  ArrayXd xs(4);
  xs << 1, 5, 6, 9;
  ArrayXd xsqr = xs * xs;
  ArrayXd ys = InterpolateUneven(xs, xsqr, ArrayXd::LinSpaced(101, 0, 10));
  float epsilon = 0.005;
  EXPECT_NEAR(1, ys(10), epsilon);
  EXPECT_NEAR(13, ys(30), epsilon);
  EXPECT_NEAR(25, ys(50), epsilon);
  EXPECT_NEAR(66, ys(80), epsilon);
}

TEST_F(FunctionOpsTest, Interpolote2DTest) {
  ArrayXd xs = ArrayXd::LinSpaced(151, 0, 10);
  ArrayXd ys = ArrayXd::LinSpaced(151, 0, 10);
  VectorXd xsqr = xs.square().matrix();
  VectorXd ysqr = ys.square().matrix();
  ArrayXXd zs = xsqr.transpose().replicate(ysqr.size(), 1) -
                ysqr.replicate(1, xsqr.size());
  ArrayXd xs_test = ArrayXd::LinSpaced(21, 0, 10);
  ArrayXd ys_test = ArrayXd::LinSpaced(21, 0, 10);
  ArrayXXd zs_test = Interpolate2D(xs, ys, zs, xs_test, ys_test);
  float epsilon = 0.005;
  EXPECT_NEAR(0, zs_test(0, 0), epsilon);
  EXPECT_NEAR(0, zs_test(10, 10), epsilon);
  EXPECT_NEAR(21, zs_test(4, 10), epsilon);
  EXPECT_NEAR(-8.75, zs_test(18, 17), epsilon);
  EXPECT_NEAR(-80, zs_test(18, 2), epsilon);
}

TEST_F(FunctionOpsTest, InterpoloteUneven2DTest) {
  ArrayXd xs(4);
  xs << 1, 5, 6, 9;
  ArrayXd ys(3);
  ys << 0, 2, 6;
  VectorXd xsqr = xs.square().matrix();
  VectorXd ysqr = ys.square().matrix();
  ArrayXXd zs = xsqr.transpose().replicate(ysqr.size(), 1) -
                ysqr.replicate(1, xsqr.size());
  ArrayXd xs_test = ArrayXd::LinSpaced(11, 0, 10);
  ArrayXd ys_test = ArrayXd::LinSpaced(11, 0, 10);
  ArrayXXd zs_test = Interpolate2DUneven(xs, ys, zs, xs_test, ys_test);
  float epsilon = 0.005;
  EXPECT_NEAR(1, zs_test(0, 0), epsilon);
  EXPECT_NEAR(0, zs_test(6, 6), epsilon);
  EXPECT_NEAR(21, zs_test(2, 5), epsilon);
  EXPECT_NEAR(45, zs_test(10, 10), epsilon);
  EXPECT_NEAR(-17, zs_test(10, 4), epsilon);
  EXPECT_NEAR(31, zs_test(4, 7), epsilon);
}


TEST_F(FunctionOpsTest, InverseUnevenTest) {
  auto [inv_xs, inv_ys] = numericaldists::InverseUneven(xs, sqr);
  ArrayXd test_xs(4);
  test_xs << 0.01, 1, 2, 80;
  ArrayXd test_ys = InterpolateUneven(inv_xs, inv_ys, test_xs);
  float epsilon = 0.001;
  EXPECT_NEAR(0.1, test_ys(0), epsilon);
  EXPECT_NEAR(1, test_ys(1), epsilon);
  EXPECT_NEAR(1.414, test_ys(2), epsilon);
  EXPECT_NEAR(8.944, test_ys(3), epsilon);
}

TEST_F(FunctionOpsTest, InverseUnevenVerticalTest) {
  float epsilon = 0.001;
  ArrayXd xs_vert(4);
  xs_vert << 20, 30, 30, 40;
  ArrayXd ys_vert(4);
  ys_vert << 5, 10, 20, 25;
  auto [inv_xs, inv_ys] = numericaldists::InverseUneven(xs_vert, ys_vert);
  ArrayXd test_xs(4);
  test_xs << 0, 15, 22, 30;
  ArrayXd test_ys = InterpolateUneven(inv_xs, inv_ys, test_xs);
  EXPECT_NEAR(20, test_ys(0), epsilon);
  EXPECT_NEAR(30, test_ys(1), epsilon);
  EXPECT_NEAR(34, test_ys(2), epsilon);
  EXPECT_NEAR(40, test_ys(3), epsilon);
}

TEST_F(FunctionOpsTest, InverseHorizontalFunctionTest) {
  float epsilon = 0.001;
  auto [inv_xs, inv_ys] = InverseUneven(xs, horz);
  ArrayXd test_xs(3);
  test_xs << 0, 5, 10;
  ArrayXd test_ys = InterpolateUneven(inv_xs, inv_ys, test_xs);
  EXPECT_NEAR(0, test_ys(0), epsilon);
  EXPECT_NEAR(10, test_ys(2), epsilon);
}

TEST_F(FunctionOpsTest, InverseDecreasingFunctionTest) {
  float epsilon = 0.001;
  ArrayXd xs = ArrayXd::LinSpaced(101, 0, 30);
  ArrayXd ys = ArrayXd::LinSpaced(101, 30, 10);
  auto [inv_xs, inv_ys] = InverseUneven(xs, ys);
  ArrayXd test_xs(3);
  test_xs << 10, 26, 30;
  ArrayXd test_ys = InterpolateUneven(inv_xs, inv_ys, test_xs);
  EXPECT_NEAR(30, test_ys(0), epsilon);
  EXPECT_NEAR(6, test_ys(1), epsilon);
  EXPECT_NEAR(0, test_ys(2), epsilon);
}

TEST_F(FunctionOpsTest, DerivativeTest) {
  float epsilon = 0.001;
  ArrayXd derivs = Derivative(xs, sqr);
  ArrayXd test_xs(3);
  test_xs << 0, 3, 8;
  ArrayXd test_ys = Interpolate(xs, derivs, test_xs);
  EXPECT_NEAR(0, test_ys(0), epsilon);
  EXPECT_NEAR(6, test_ys(1), epsilon);
  EXPECT_NEAR(16, test_ys(2), epsilon);
}

TEST_F(FunctionOpsTest, AreasTest) {
  float epsilon = 0.001;
  ArrayXd xs(4);
  xs << 1, 4, 5, 9;
  ArrayXd areas = Areas(xs, xs * xs);
  EXPECT_NEAR(0, areas(0), epsilon);
  EXPECT_NEAR(51.0/2, areas(1), epsilon);
  EXPECT_NEAR(41.0/2, areas(2), epsilon);
  EXPECT_NEAR(212.0, areas(3), epsilon);  
}

TEST_F(FunctionOpsTest, Areas2DTest) {
  ArrayXd xs(4);
  xs << 1, 5, 6, 9;
  ArrayXd ys(3);
  ys << 0, 2, 5;
  VectorXd xsqr = xs.square().matrix();
  VectorXd ysqr = ys.square().matrix();
  ArrayXXd zs = xsqr.transpose().replicate(ysqr.size(), 1) -
                ysqr.replicate(1, xsqr.size());
  ArrayXXd areas = Areas2D(xs, ys, zs);
  float epsilon = 0.001;
  EXPECT_NEAR(0, areas(0, 0), epsilon);
  EXPECT_NEAR(0, areas(2, 0), epsilon);
  EXPECT_NEAR(0, areas(0, 3), epsilon);
  EXPECT_NEAR(9*(234-58)/4.0, areas(2, 3), epsilon);  
}

TEST_F(FunctionOpsTest, IntegralBelowTest) {
  float epsilon = 0.001;
  ArrayXd xs = ArrayXd::LinSpaced(101, -1, 1);
  ArrayXd integral = IntegralBelow(xs, xs * xs);
  ArrayXd test_xs(3);
  test_xs << -1, 0, 1;
  ArrayXd test_ys = Interpolate(xs, integral, test_xs);
  EXPECT_NEAR(0, test_ys(0), epsilon);
  EXPECT_NEAR(0.3333, test_ys(1), epsilon);
  EXPECT_NEAR(0.6666, test_ys(2), epsilon);
}

TEST_F(FunctionOpsTest, IntegralBelow2DTest) {
  ArrayXd xs(4);
  xs << 1, 5, 6, 9;
  ArrayXd ys(3);
  ys << 0, 2, 5;
  VectorXd xsqr = xs.square().matrix();
  VectorXd ysqr = ys.square().matrix();
  ArrayXXd zs = xsqr.transpose().replicate(ysqr.size(), 1) -
                ysqr.replicate(1, xsqr.size());
  ArrayXXd areas = IntegralBelow2D(xs, ys, zs);
  float epsilon = 0.001;
  EXPECT_NEAR(0, areas(0, 0), epsilon);
  EXPECT_NEAR(0, areas(2, 0), epsilon);
  EXPECT_NEAR(0, areas(0, 3), epsilon);
  EXPECT_NEAR(88+57, areas(1, 2), epsilon);
  EXPECT_NEAR(88+57+30, areas(2, 2), epsilon);
}

TEST_F(FunctionOpsTest, IntegralAboveTest) {
  float epsilon = 0.001;
  ArrayXd xs = ArrayXd::LinSpaced(101, -1, 1);
  ArrayXd integral = IntegralAbove(xs, xs * xs);
  ArrayXd test_xs(3);
  test_xs << -1, 0, 1;
  ArrayXd test_ys = Interpolate(xs, integral, test_xs);
  EXPECT_NEAR(0.6666, test_ys(0), epsilon);
  EXPECT_NEAR(0.3333, test_ys(1), epsilon);
  EXPECT_NEAR(0, test_ys(2), epsilon);
}

}  // namespace gatests
