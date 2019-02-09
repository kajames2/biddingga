#include <gtest/gtest.h>
#include <vector>
#include "auctions/common_value_signal.h"
#include "boost/math/distributions/uniform.hpp"
#include "numericaldists/distribution.h"
#include "numericaldists/function_ops.h"
#include "numericaldists/piecewise_linear.h"
namespace gatests {

class CommonValueSignalTest : public ::testing::Test {
 public:
  CommonValueSignalTest() {}

 protected:
  virtual void SetUp() {}
  auctions::CommonValueSignal auction = auctions::CommonValueSignal({2, 2}, 1, {-4, 4});
  auctions::CommonValueSignal auction_single = auctions::CommonValueSignal({1, 1}, 1, {-4, 4});
  auctions::CommonValueSignal auction_mixed = auctions::CommonValueSignal({1, 2}, 1, {-4, 4});
  float epsilon = 0.0001;
};

TEST_F(CommonValueSignalTest, AlwaysWinTest) {
  auction.AcceptStrategy(numericaldists::PiecewiseLinear({3, 3},{0, 2}), 0);
  auction.AcceptStrategy(numericaldists::PiecewiseLinear({1, 1},{0, 2}), 1);
  float fit = auction.GetFitness(numericaldists::PiecewiseLinear({1, 1},{0, 2}), 1);
  EXPECT_NEAR(1, fit, epsilon);
}

TEST_F(CommonValueSignalTest, WinHalfTest) {
  auto bid_func = numericaldists::PiecewiseLinear({0, 0},{0, 2});
  auction.AcceptStrategy(bid_func, 0);
  auction.AcceptStrategy(bid_func, 1);
  float fit = auction.GetFitness(bid_func, 1);
  EXPECT_NEAR(-0.116666, fit, epsilon);
}


TEST_F(CommonValueSignalTest, AlwaysWinSinglesTest) {
  auction_single.AcceptStrategy(3, 0);
  auction_single.AcceptStrategy(1, 1);
  float fit = auction_single.GetFitness(1, 1);
  EXPECT_NEAR(1, fit, epsilon);
}

TEST_F(CommonValueSignalTest, WinHalfSinglesTest) {
  auction_single.AcceptStrategy(0, 0);
  auction_single.AcceptStrategy(0, 1);
  float fit = auction_single.GetFitness(0, 1);
  EXPECT_NEAR(-0.166666, fit, epsilon);
}


TEST_F(CommonValueSignalTest, AlwaysWinMixedSingleTest) {
  auction_mixed.AcceptStrategy(3, 0);
  auction_mixed.AcceptStrategy(numericaldists::PiecewiseLinear({1, 1},{0, 2}), 1);
  float fit = auction_mixed.GetFitness(numericaldists::PiecewiseLinear({1, 1},{0, 2}), 1);
  EXPECT_NEAR(1, fit, epsilon);
}

TEST_F(CommonValueSignalTest, AlwaysWinMixedMultiTest) {
  auction_mixed.AcceptStrategy(1, 0);
  auction_mixed.AcceptStrategy(numericaldists::PiecewiseLinear({3, 3},{0, 2}), 1);
  float fit = auction_mixed.GetFitness(1, 0);
  EXPECT_NEAR(1, fit, epsilon);
}


}  // namespace gatests
