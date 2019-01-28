#include <gtest/gtest.h>

#include <vector>

#include "numericaldists/uneven_piecewise_linear.h"

namespace gatests {

using namespace numericaldists;

class UnevenPiecewiseLinearTest : public ::testing::Test {
 public:
  UnevenPiecewiseLinearTest() {}

 protected:
  virtual void SetUp() {}
  UnevenPiecewiseLinear func = UnevenPiecewiseLinear(
      std::vector<float>{0, 10, 15, 30}, std::vector<float>{10, 20, 20, 50});
  UnevenPiecewiseLinear horizontal_func = UnevenPiecewiseLinear(
      std::vector<float>{0, 30}, std::vector<float>{10, 10});
  UnevenPiecewiseLinear vertical_func = UnevenPiecewiseLinear(
      std::vector<float>{30, 30}, std::vector<float>{10, 20});
  UnevenPiecewiseLinear decreasing_func = UnevenPiecewiseLinear(
      std::vector<float>{0, 30}, std::vector<float>{30, 10});
};

typedef UnevenPiecewiseLinearTest UnevenPiecewiseLinearDeathTest;

TEST_F(UnevenPiecewiseLinearTest, GetBidInterior) {
  EXPECT_EQ(15, func.GetBid(5));
  EXPECT_EQ(20, func.GetBid(10));
  EXPECT_EQ(20, func.GetBid(13));
  EXPECT_EQ(40, func.GetBid(25));
}

TEST_F(UnevenPiecewiseLinearTest, GetBidExtrapolation) {
  EXPECT_EQ(10, func.GetBid(0));
  EXPECT_EQ(50, func.GetBid(30));
  EXPECT_EQ(50, func.GetBid(50));
  EXPECT_EQ(10, func.GetBid(-120));
}

TEST_F(UnevenPiecewiseLinearTest, GetBidVertical) {
  EXPECT_EQ(10, vertical_func.GetBid(0));
  EXPECT_EQ(20, vertical_func.GetBid(50));
}

TEST_F(UnevenPiecewiseLinearTest, GetBidDecreasing) {
  EXPECT_FLOAT_EQ(30, decreasing_func.GetBid(-10));
  EXPECT_FLOAT_EQ(26, decreasing_func.GetBid(6));
  EXPECT_FLOAT_EQ(10, decreasing_func.GetBid(40));
}

TEST_F(UnevenPiecewiseLinearDeathTest, SinglePointFunctionDeathTest) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH(UnevenPiecewiseLinear(std::vector<float>(0, 30), std::vector<float>{10}),
               "");
}

}  // namespace gatests
