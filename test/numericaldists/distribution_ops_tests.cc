#include <gtest/gtest.h>

#include <vector>

#include <boost/math/distributions/normal.hpp>
#include "numericaldists/distribution.h"
#include "numericaldists/distribution_ops.h"

namespace gatests {

using namespace numericaldists;

auto lin = [](float x) { return 2 * x; };
auto sqr = [](float x) { return 1.5 * x * x; };
auto cub = [](float x) { return 0.5 * (1 + x * x * x); };

class DistributionOpsTest : public ::testing::Test {
 public:
  DistributionOpsTest() {}

 protected:
  virtual void SetUp() {}
};

TEST_F(DistributionOpsTest, ApproximatePDFExpectedValueFunctionTest) {
  auto func = ApproximatePDFExpectedValueFunction(sqr, Interval{-1, 1});
  float epsilon = 0.001;
  EXPECT_NEAR(0, func(1), epsilon);
  EXPECT_NEAR(-3 / 4., func(0), epsilon);
  EXPECT_NEAR(-1, func(-1), epsilon);
}

TEST_F(DistributionOpsTest, ApproximateCDFExpectedValueFunctionTest) {
  auto func = ApproximateCDFExpectedValueFunction(cub, Interval{-1, 1});
  float epsilon = 0.001;
  EXPECT_NEAR(0, func(1), epsilon);
  EXPECT_NEAR(-3 / 4., func(0), epsilon);
  EXPECT_NEAR(-1, func(-1), epsilon);
}

TEST_F(DistributionOpsTest, ApproximateRandomVariableFunctionCDFTest) {
  Distribution dist(boost::math::normal_distribution<>(0, 1));
  auto func = ApproximateRandomVariableFunctionCDF(dist, lin);
  float epsilon = 0.001;
  EXPECT_NEAR(1 / 2., func(0), epsilon);
  EXPECT_NEAR(0.841345, func(2), epsilon);
  EXPECT_NEAR(0.022750, func(-4), epsilon);
}

TEST_F(DistributionOpsTest, ApproximateRandomVariableFunctionCDF2DTest) {
  auto pdf = [](float x, float y) { return x / 2 + 2 * y; };
  auto func = [](float x, float y) { return 1 / 4. * x * x - 4 * y * y; };
  auto out = ApproximateRandomVariableFunctionCDF(pdf, func, {0, 2}, {0, 0.5},
                                                  {-1, 1});
  float epsilon = 0.001;
  EXPECT_NEAR(0, out(-2), epsilon);
  EXPECT_NEAR(0, out(-1), epsilon);
  EXPECT_NEAR(0.152369, out(-.5), epsilon);
  EXPECT_NEAR(1 / 2., out(0), epsilon);
  EXPECT_NEAR(1, out(1), epsilon);
  EXPECT_NEAR(1, out(2), epsilon); 
}

// std::function<float(float, float)> MultivariatePDFDomainTransform(
//     const std::function<float(float, float)>& func,
//     const std::function<float(float, float)>& x_trans,
//     const std::function<float(float, float)>& y_trans,
//     const std::function<float(float, float)>& jacobian);

}  // namespace gatests
