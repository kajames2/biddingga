#include <gtest/gtest.h>

#include <vector>

#include "numericaldists/function_ops.h"
#include "numericaldists/piecewise_linear.h"

namespace gatests {

using namespace numericaldists;

auto horz = [] (float x) { return 10; };
auto lin = [] (float x) { return 2 * x; };
auto sqr = [] (float x) { return x * x; };
auto cube = [] (float x) { return x * x * x; };

class FunctionOpsTest : public ::testing::Test {
 public:
  FunctionOpsTest() {}
 protected:
  virtual void SetUp() {}
};

TEST_F(FunctionOpsTest, SampleFunctionTest) {
  auto func =
      ResampleFunction(sqr, Interval{0, 10}, 101);
  float epsilon = 0.005;
  EXPECT_NEAR(0.01, func(0.1), epsilon);
  EXPECT_NEAR(1, func(1), epsilon);
  EXPECT_NEAR(2, func(1.414), epsilon);
  EXPECT_NEAR(80, func(8.944), epsilon);
}

TEST_F(FunctionOpsTest, ApproximateInverseTest) {
  auto inv = ApproximateInverse(sqr, Interval{0, 10});
  float epsilon = 0.001;
  EXPECT_NEAR(0.1, inv(0.01), epsilon);
  EXPECT_NEAR(1, inv(1), epsilon);
  EXPECT_NEAR(1.414, inv(2), epsilon);
  EXPECT_NEAR(8.944, inv(80), epsilon);
}

TEST_F(FunctionOpsTest, InverseVerticalFunctionTest) {
  float epsilon = 0.001;
  auto vertical_func =
      PiecewiseLinear(std::vector<float>{10, 20}, Interval{30, 30});
  auto func_inv = ApproximateInverse(vertical_func, Interval{30, 30});
  EXPECT_NEAR(30, func_inv(-10), epsilon);
  EXPECT_NEAR(30, func_inv(15), epsilon);
  EXPECT_NEAR(30, func_inv(30), epsilon);
}

TEST_F(FunctionOpsTest, InverseHorizontalFunctionTest) {
  float epsilon = 0.001;
  auto func_inv = ApproximateInverse(horz, Interval{0, 30});
  EXPECT_NEAR(0, func_inv(0), epsilon);
  EXPECT_NEAR(30, func_inv(20), epsilon);
}

TEST_F(FunctionOpsTest, InverseDecreasingFunctionTest) {
  float epsilon = 0.001;
  auto decreasing_func =
      PiecewiseLinear(std::vector<float>{30, 10}, Interval{0, 30});
  auto func_inv = ApproximateInverse(decreasing_func, Interval{0, 30});
  EXPECT_NEAR(0, func_inv(30), epsilon);
  EXPECT_NEAR(6, func_inv(26), epsilon);
  EXPECT_NEAR(30, func_inv(10), epsilon);
}

TEST_F(FunctionOpsTest, ApproximateDerivativeTest) {
  float epsilon = 0.001;
  auto deriv = ApproximateDerivative(sqr, Interval{-5, 5});
  EXPECT_NEAR(0, deriv(0), epsilon);
  EXPECT_NEAR(-6, deriv(-3), epsilon);
  EXPECT_NEAR(6, deriv(3), epsilon);
}

TEST_F(FunctionOpsTest, ApproximateIntegralBelowTest) {
  float epsilon = 0.001;
  auto integral = ApproximateIntegralBelow(sqr, Interval{-1,1});
  EXPECT_NEAR(0, integral(-1), epsilon);
  EXPECT_NEAR(0.3333, integral(0), epsilon);
  EXPECT_NEAR(0.6666, integral(1), epsilon);
}

TEST_F(FunctionOpsTest, ApproximateIntegralAboveTest) {
  float epsilon = 0.001;
  auto integral = ApproximateIntegralAbove(sqr, Interval{-1,1});
  EXPECT_NEAR(0, integral(1), epsilon);
  EXPECT_NEAR(0.3333, integral(0), epsilon);
  EXPECT_NEAR(0.6666, integral(-1), epsilon);
}

TEST_F(FunctionOpsTest, SampleFunction2DTest) {
  auto func =
      ResampleFunction2D([](float x, float y) {return x * x - y * y;}, Interval{0, 10}, Interval{0,10});
  float epsilon = 0.005;
  EXPECT_NEAR(0, func(0, 0), epsilon);
  EXPECT_NEAR(0, func(5, 5), epsilon);
  EXPECT_NEAR(21, func(5 ,2), epsilon);
  EXPECT_NEAR(-8.75, func(8.5, 9), epsilon);
  EXPECT_NEAR(-80, func(1, 9), epsilon);
}

}  // namespace gatests
