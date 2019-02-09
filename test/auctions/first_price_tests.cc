#include <gtest/gtest.h>
#include <vector>
#include "numericaldists/function_ops.h"
#include "numericaldists/distribution.h"
#include "auctions/first_price.h"
#include "boost/math/distributions/uniform.hpp"

namespace gatests {

class FirstPriceTest : public ::testing::Test {
 public:
  FirstPriceTest() {}

 protected:
  virtual void SetUp() {}
  auctions::FirstPrice auction =
      auctions::FirstPrice(std::vector<numericaldists::Distribution>{
          boost::math::uniform_distribution<>(0, 1),
          boost::math::uniform_distribution<>(0, 1)});
  float epsilon = 0.0001;
};

TEST_F(FirstPriceTest, AlwaysWinTest) {
  auction.AcceptStrategy([](float x) { return 0; }, 0);
  float fit = auction.GetFitness([](float x) { return 0.01; }, 1);
  EXPECT_NEAR(0.49, fit, epsilon);
}

TEST_F(FirstPriceTest, AlwaysTieTest) {
  auction.AcceptStrategy([](float x) { return 0; }, 0);
  float fit = auction.GetFitness([](float x) { return 0; }, 1);
  EXPECT_NEAR(0.25, fit, epsilon);
}

TEST_F(FirstPriceTest, WinHalfTest) {
  auto bid_func = [](float x) { return x / 2; };
  auction.AcceptStrategy(bid_func, 0);
  float fit = auction.GetFitness([](float x) { return x / 2; }, 1);
  EXPECT_NEAR(.166666, fit, epsilon);
}


TEST_F(FirstPriceTest, DecreasingBidsTest) {
  auto bid_func = [](float x) { return 0.5 - x / 2; };
  auction.AcceptStrategy(bid_func, 0);
  float fit = auction.GetFitness([](float x) { return x / 2; }, 1);
  EXPECT_NEAR(.166666, fit, epsilon);
}

TEST_F(FirstPriceTest, DecreasingBids2Test) {
  auto bid_func = [](float x) { return x / 2; };
  auction.AcceptStrategy(bid_func, 0);
  float fit = auction.GetFitness([](float x) { return 0.5 - 0.5*x; }, 1);
  EXPECT_NEAR(0, fit, epsilon);
}

TEST_F(FirstPriceTest, NoIntersectDistsTest) {
  auction = auctions::FirstPrice(std::vector<numericaldists::Distribution>{
      boost::math::uniform_distribution<>(20, 40),
      boost::math::uniform_distribution<>(40, 60)});
  float epsilon = 0.001;
  auto bid_func = [](float x) { return 20; };
  auction.AcceptStrategy(bid_func, 0);
  auto bid_func2 = [](float x) { return 20.2; };
  auction.AcceptStrategy(bid_func2, 1);
  float fit = auction.GetFitness(bid_func, 0);
  EXPECT_NEAR(0, fit, epsilon);
  float fit2 = auction.GetFitness(bid_func2, 1);
  EXPECT_NEAR(29.8, fit2, epsilon);
}

}  // namespace gatests
