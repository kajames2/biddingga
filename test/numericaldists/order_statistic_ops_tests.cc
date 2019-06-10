#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include <boost/math/distributions/normal.hpp>
#include "numericaldists/function_ops.h"
#include "numericaldists/order_statistic_ops.h"

#include <eigen3/Eigen/Core>

namespace gatests {

using namespace Eigen;
using namespace numericaldists;

int n_samples = 3001;

ArrayXd xs2 = ArrayXd::LinSpaced(n_samples, -1, 2);
auto uni_cdf = [](double x) -> double {
  if (x < 0) {
    return 0;
  } else if (x > 2) {
    return 1;
  } else {
    return x / 2;
  }
};

ArrayXd uni_cdf_s = xs2.unaryExpr(uni_cdf);

auto lin_cdf = [](double x) -> double {
  if (x < 0) {
    return 0;
  } else if (x > 1) {
    return 1;
  } else {
    return x * x;
  }
};
ArrayXd lin_cdf_s = xs2.unaryExpr(lin_cdf);

auto sqr_cdf = [](double x) -> double {
  if (x < -1) {
    return 0;
  } else if (x > 1) {
    return 1;
  } else {
    return 0.5 * (1 + x * x * x);
  }
};
ArrayXd sqr_cdf_s = xs2.unaryExpr(sqr_cdf);

auto uni_cdf_ord_3_1 = [](double x) {
  return std::pow(x / 2, 3) - 3 * std::pow(x / 2, 2) + 3 * x / 2;
};
auto uni_cdf_ord_3_2 = [](double x) {
  return 3 * std::pow(x / 2, 2) - 2 * std::pow(x / 2, 3);
};
auto uni_cdf_ord_3_3 = [](double x) { return std::pow(x / 2, 3); };

auto multi_cdf_ord_1_neg = [](double x) { return 0.5 * (1 + x * x * x); };
auto multi_cdf_ord_1_pos = [](double x) {
  return 0.25 * (2 + x + 2 * std::pow(x, 2) + std::pow(x, 3) - std::pow(x, 4) -
                 2 * std::pow(x, 5) + std::pow(x, 6));
};
auto multi_cdf_ord_2 = [](double x) {
  return 0.25 * (x + 2 * std::pow(x, 2) + std::pow(x, 4) + 2 * std::pow(x, 5) -
                 2 * std::pow(x, 6));
};
auto multi_cdf_ord_3_01 = [](double x) {
  return 0.25 * std::pow(x, 3) * (1 + std::pow(x, 3));
};
auto multi_cdf_ord_3_12 = [](double x) { return 0.5 * x; };

class OrderStatisticOpsTest : public ::testing::Test {
 public:
  OrderStatisticOpsTest() {}

 protected:
  virtual void SetUp() {}
};

TEST_F(OrderStatisticOpsTest, LowestSingleDistOrderStatisticTest) {
  ArrayXd order_stats = KthLowestOrderStatisticCDF(xs2, uni_cdf_s, 3, 1);
  ArrayXd test_xs(3);
  test_xs << 0, 0.75, 2;
  ArrayXd test_ys = Interpolate(xs2, order_stats, test_xs);

  float epsilon = 0.001;
  EXPECT_NEAR(0, test_ys(0), epsilon);
  EXPECT_NEAR(uni_cdf_ord_3_1(0.75), test_ys(1), epsilon);
  EXPECT_NEAR(1, test_ys(2), epsilon);
}

TEST_F(OrderStatisticOpsTest, MiddleSingleDistOrderStatisticTest) {
  ArrayXd order_stats = KthLowestOrderStatisticCDF(xs2, uni_cdf_s, 3, 2);
  ArrayXd test_xs(3);
  test_xs << 0, 0.75, 2;
  ArrayXd test_ys = Interpolate(xs2, order_stats, test_xs);

  float epsilon = 0.001;
  EXPECT_NEAR(0, test_ys(0), epsilon);
  EXPECT_NEAR(uni_cdf_ord_3_2(0.75), test_ys(1), epsilon);
  EXPECT_NEAR(1, test_ys(2), epsilon);
}

TEST_F(OrderStatisticOpsTest, HighestSingleDistOrderStatisticTest) {
  ArrayXd order_stats = KthLowestOrderStatisticCDF(xs2, uni_cdf_s, 3, 3);
  ArrayXd test_xs(3);
  test_xs << 0, 0.75, 2;
  ArrayXd test_ys = Interpolate(xs2, order_stats, test_xs);

  float epsilon = 0.001;
  EXPECT_NEAR(0, test_ys(0), epsilon);
  EXPECT_NEAR(uni_cdf_ord_3_3(0.75), test_ys(1), epsilon);
  EXPECT_NEAR(1, test_ys(2), epsilon);
}

TEST_F(OrderStatisticOpsTest, LowestMultiDistOrderStatisticTest) {
  ArrayXd order_stats =
      KthLowestOrderStatisticCDF(xs2, {uni_cdf_s, lin_cdf_s, sqr_cdf_s}, 1);
  ArrayXd test_xs(5);
  test_xs << -1, -0.75, 0.75, 1, 1.5;
  ArrayXd test_ys = Interpolate(xs2, order_stats, test_xs);

  float epsilon = 0.001;
  EXPECT_NEAR(0, test_ys(0), epsilon);
  EXPECT_NEAR(multi_cdf_ord_1_neg(-0.75), test_ys(1), epsilon);
  EXPECT_NEAR(multi_cdf_ord_1_pos(0.75), test_ys(2), epsilon);
  EXPECT_NEAR(1, test_ys(3), epsilon);
  EXPECT_NEAR(1, test_ys(4), epsilon);
}

TEST_F(OrderStatisticOpsTest, MiddleMultiDistOrderStatisticTest) {
  ArrayXd order_stats =
      KthLowestOrderStatisticCDF(xs2, {uni_cdf_s, lin_cdf_s, sqr_cdf_s}, 2);
  ArrayXd test_xs(5);
  test_xs << -1, 0, 0.75, 1, 1.5;
  ArrayXd test_ys = Interpolate(xs2, order_stats, test_xs);

  float epsilon = 0.001;
  EXPECT_NEAR(0, test_ys(0), epsilon);
  EXPECT_NEAR(0, test_ys(1), epsilon);
  EXPECT_NEAR(multi_cdf_ord_2(0.75), test_ys(2), epsilon);
  EXPECT_NEAR(1, test_ys(3), epsilon);
  EXPECT_NEAR(1, test_ys(4), epsilon);
}

// TEST_F(OrderStatisticOpsTest, MiddleMultiDistOrderStatisticTest) {
//   auto func =
//       KthLowestOrderStatisticCDF({uni_cdf_s, lin_cdf_s, sqr_cdf_s}, {-1, 2},
//       2);
//   float epsilon = 0.001;
//   EXPECT_NEAR(1, func(1.5), epsilon);
//   EXPECT_NEAR(1, func(1), epsilon);
//   EXPECT_NEAR(multi_cdf_ord_2(0.75), func(0.75), epsilon);
//   EXPECT_NEAR(0, func(0), epsilon);
//   EXPECT_NEAR(0, func(-1), epsilon);
// }

TEST_F(OrderStatisticOpsTest, HighestMultiDistOrderStatisticTest) {
  ArrayXd order_stats =
      KthLowestOrderStatisticCDF(xs2, {uni_cdf_s, lin_cdf_s, sqr_cdf_s}, 3);
  ArrayXd test_xs(6);
  test_xs << -0.5, 0, 0.75, 1, 1.5, 2.5;
  ArrayXd test_ys = Interpolate(xs2, order_stats, test_xs);

  float epsilon = 0.001;
  EXPECT_NEAR(0, test_ys(0), epsilon);
  EXPECT_NEAR(0, test_ys(1), epsilon);
  EXPECT_NEAR(multi_cdf_ord_3_01(0.75), test_ys(2), epsilon);
  EXPECT_NEAR(multi_cdf_ord_3_01(1), test_ys(3), epsilon);
  EXPECT_NEAR(multi_cdf_ord_3_12(1.5), test_ys(4), epsilon);
  EXPECT_NEAR(1, test_ys(5), epsilon);
}


// TEST_F(OrderStatisticOpsTest, HighestMultiDistOrderStatisticTest) {
//   auto func =
//       KthLowestOrderStatisticCDF({uni_cdf_s, lin_cdf_s, sqr_cdf_s}, {-1, 2},
//       3);
//   float epsilon = 0.001;
//   EXPECT_NEAR(1, func(2.5), epsilon);
//   EXPECT_NEAR(multi_cdf_ord_3_12(1.5), func(1.5), epsilon);
//   EXPECT_NEAR(multi_cdf_ord_3_01(1), func(1), epsilon);
//   EXPECT_NEAR(multi_cdf_ord_3_01(0.75), func(0.75), epsilon);
//   EXPECT_NEAR(0, func(0), epsilon);
//   EXPECT_NEAR(0, func(-0.5), epsilon);
// }

}  // namespace gatests
