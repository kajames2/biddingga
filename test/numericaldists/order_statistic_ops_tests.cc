#include <gtest/gtest.h>

#include <vector>
#include <cmath>

#include <boost/math/distributions/normal.hpp>
#include "numericaldists/distribution.h"
#include "numericaldists/function_ops.h"
#include "numericaldists/order_statistic_ops.h"

namespace gatests {

using namespace numericaldists;

auto uni_cdf = [](float x) { return x/2; };
Interval uni_int{0, 2};
auto uni_cdf_s = ResampleFunction(uni_cdf, uni_int);

auto lin_cdf = [](float x) { return x * x; };
Interval lin_int{0, 1};
auto lin_cdf_s = ResampleFunction(lin_cdf, lin_int);

auto sqr_cdf = [](float x) { return 0.5 * (1 + x * x * x); };
Interval sqr_int{-1, 1};
auto sqr_cdf_s = ResampleFunction(sqr_cdf, sqr_int);

auto uni_cdf_ord_3_1 = [](float x) { return std::pow(x/2,3) - 3 * std::pow(x/2,2) + 3 * x/2; };
auto uni_cdf_ord_3_2 = [](float x) { return 3 * std::pow(x/2,2) - 2 * std::pow(x/2,3); };
auto uni_cdf_ord_3_3 = [](float x) { return std::pow(x/2,3); };

auto multi_cdf_ord_1_neg = [](float x) {return 0.5 * (1 + x*x*x);} ;
auto multi_cdf_ord_1_pos = [](float x) {return 0.25 * (2 + x + 2 * std::pow(x,2) + std::pow(x,3) - std::pow(x,4) - 2 * std::pow(x,5) + std::pow(x,6));};
auto multi_cdf_ord_2 = [](float x) {return 0.25 * (x + 2*std::pow(x,2) + std::pow(x,4) + 2*std::pow(x,5) - 2 * std::pow(x,6));};
auto multi_cdf_ord_3_01 = [](float x) {return 0.25 * std::pow(x,3) * (1 + std::pow(x,3));};
auto multi_cdf_ord_3_12 = [](float x) {return 0.5*x;};


class OrderStatisticOpsTest : public ::testing::Test {
 public:
  OrderStatisticOpsTest() {}

 protected:
  virtual void SetUp() {}
};

TEST_F(OrderStatisticOpsTest, LowestSingleDistOrderStatisticTest) {
  auto func = ApproximateKthLowestOrderStatisticCDF(uni_cdf, uni_int, 3, 1);
  float epsilon = 0.001;
  EXPECT_NEAR(1, func(2), epsilon);
  EXPECT_NEAR(uni_cdf_ord_3_1(0.75), func(0.75), epsilon);
  EXPECT_NEAR(0, func(0), epsilon);
}

TEST_F(OrderStatisticOpsTest, MiddleSingleDistOrderStatisticTest) {
  auto func = ApproximateKthLowestOrderStatisticCDF(uni_cdf, uni_int, 3, 2);
  float epsilon = 0.001;
  EXPECT_NEAR(1, func(2), epsilon);
  EXPECT_NEAR(uni_cdf_ord_3_2(0.75), func(0.75), epsilon);
  EXPECT_NEAR(0, func(0), epsilon);
}

TEST_F(OrderStatisticOpsTest, HighestSingleDistOrderStatisticTest) {
  auto func = ApproximateKthLowestOrderStatisticCDF(uni_cdf, uni_int, 3, 3);
  float epsilon = 0.001;
  EXPECT_NEAR(1, func(2), epsilon);
  EXPECT_NEAR(uni_cdf_ord_3_3(0.75), func(0.75), epsilon);
  EXPECT_NEAR(0, func(0), epsilon);
}

TEST_F(OrderStatisticOpsTest, LowestMultiDistOrderStatisticTest) {
  auto func = ApproximateKthLowestOrderStatisticCDF({uni_cdf_s, lin_cdf_s, sqr_cdf_s},
                                                    {-1,2}, 1);
  float epsilon = 0.001;
  EXPECT_NEAR(1, func(1.5), epsilon);
  EXPECT_NEAR(1, func(1), epsilon);
  EXPECT_NEAR(multi_cdf_ord_1_pos(0.75), func(0.75), epsilon);
  EXPECT_NEAR(multi_cdf_ord_1_neg(-0.75), func(-0.75), epsilon);
  EXPECT_NEAR(0, func(-1), epsilon);
}

TEST_F(OrderStatisticOpsTest, MiddleMultiDistOrderStatisticTest) {
  auto func = ApproximateKthLowestOrderStatisticCDF({uni_cdf_s, lin_cdf_s, sqr_cdf_s},
                                                    {-1,2}, 2);
  float epsilon = 0.001;
  EXPECT_NEAR(1, func(1.5), epsilon);
  EXPECT_NEAR(1, func(1), epsilon);
  EXPECT_NEAR(multi_cdf_ord_2(0.75), func(0.75), epsilon);
  EXPECT_NEAR(0, func(0), epsilon);
  EXPECT_NEAR(0, func(-1), epsilon);
}

TEST_F(OrderStatisticOpsTest, HighestMultiDistOrderStatisticTest) {
  auto func = ApproximateKthLowestOrderStatisticCDF({uni_cdf_s, lin_cdf_s, sqr_cdf_s},
                                                    {-1,2}, 3);
  float epsilon = 0.001;
  EXPECT_NEAR(1, func(2.5), epsilon);
  EXPECT_NEAR(multi_cdf_ord_3_12(1.5), func(1.5), epsilon);
  EXPECT_NEAR(multi_cdf_ord_3_01(1), func(1), epsilon);
  EXPECT_NEAR(multi_cdf_ord_3_01(0.75), func(0.75), epsilon);
  EXPECT_NEAR(0, func(0), epsilon);
  EXPECT_NEAR(0, func(-0.5), epsilon);
}

}  // namespace gatests
