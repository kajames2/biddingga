#include <gtest/gtest.h>

#include <vector>

#include "numericaldists/piecewise_linear.h"

namespace gatests {

using namespace numericaldists;

class PiecewiseLinearTest : public ::testing::Test {
 public:
  PiecewiseLinearTest() {}

 protected:
  virtual void SetUp() {}
  PiecewiseLinear func =
      PiecewiseLinear(std::vector<float>{10, 20, 20, 50}, Interval{0, 30});
  PiecewiseLinear horizontal_func =
      PiecewiseLinear(std::vector<float>{10, 10}, Interval{0, 30});
  PiecewiseLinear vertical_func =
      PiecewiseLinear(std::vector<float>{10, 20}, Interval{30, 30});
  PiecewiseLinear decreasing_func =
      PiecewiseLinear(std::vector<float>{30, 10}, Interval{0, 30});
};

typedef PiecewiseLinearTest PiecewiseLinearDeathTest;

TEST_F(PiecewiseLinearTest, GetBidInterior) {
  EXPECT_EQ(15, func.GetBid(5));
  EXPECT_EQ(20, func.GetBid(10));
  EXPECT_EQ(20, func.GetBid(15));
  EXPECT_EQ(35, func.GetBid(25));
}

TEST_F(PiecewiseLinearTest, GetBidExtrapolation) {
  EXPECT_EQ(10, func.GetBid(0));
  EXPECT_EQ(50, func.GetBid(30));
  EXPECT_EQ(50, func.GetBid(50));
  EXPECT_EQ(10, func.GetBid(-120));
}

TEST_F(PiecewiseLinearTest, GetBidVertical) {
  EXPECT_EQ(10, vertical_func.GetBid(0));
  EXPECT_EQ(20, vertical_func.GetBid(50));
}

TEST_F(PiecewiseLinearTest, GetBidDecreasing) {
  EXPECT_FLOAT_EQ(30, decreasing_func.GetBid(-10));
  EXPECT_FLOAT_EQ(26, decreasing_func.GetBid(6));
  EXPECT_FLOAT_EQ(10, decreasing_func.GetBid(40));
}

TEST_F(PiecewiseLinearDeathTest, SinglePointFunctionDeathTest) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH(PiecewiseLinear(std::vector<float>{10}, Interval(0, 30)), "");
}

}  // namespace gatests
