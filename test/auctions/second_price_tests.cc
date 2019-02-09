#include <gtest/gtest.h>
#include <vector>
#include "auctions/second_price.h"
#include "boost/math/distributions/uniform.hpp"
#include "numericaldists/function_ops.h"
#include "numericaldists/distribution.h"

namespace gatests {

class SecondPriceTest : public ::testing::Test {
 public:
  SecondPriceTest() {}

 protected:
  virtual void SetUp() {}
  auctions::SecondPrice auction =
      auctions::SecondPrice(std::vector<numericaldists::Distribution>{
          boost::math::uniform_distribution<>(0, 1),
          boost::math::uniform_distribution<>(0, 1)});
  float epsilon = 0.001;
};

// TEST_F(SecondPriceTest, AlwaysWinTest) {
//   auto func1 = [](float x) { return 0; };
//   auto func2 = [](float x) { return 1; };
//   auction.AcceptStrategy(func1, 0);
//   auction.AcceptStrategy(func2, 1);
//   float fit = auction.GetFitness(func2, 1);
//   EXPECT_NEAR(0.5, fit, epsilon);
// }

// TEST_F(SecondPriceTest, WinHalfTest) {
//   auto bid_func = [](float x) { return x; };
//   auction.AcceptStrategy(bid_func, 0);
//   auction.AcceptStrategy(bid_func, 1);
//   float fit = auction.GetFitness([](float x) { return x; }, 1);
//   EXPECT_NEAR(.166666, fit, epsilon);
// }

// TEST_F(SecondPriceTest, DecreasingBidsTest) {
//   auto bid_func = [](float x) { return 1 - x; };
//   auto bid_func2 = [](float x) { return x; };
//   auction.AcceptStrategy(bid_func, 0);
//   auction.AcceptStrategy(bid_func2, 1);
//   float fit = auction.GetFitness(bid_func2, 1);
//   EXPECT_NEAR(.166666, fit, epsilon);
// }

// TEST_F(SecondPriceTest, DecreasingBids2Test) {
//   auto bid_func = [](float x) { return x; };
//   auto bid_func2 = [](float x) { return 1 - x; };
//   auction.AcceptStrategy(bid_func, 0);
//   auction.AcceptStrategy(bid_func2, 1);
//   float fit = auction.GetFitness(bid_func2, 1);
//   EXPECT_NEAR(0, fit, epsilon);
// }

TEST_F(SecondPriceTest, NoIntersectDistsTest) {
  auction = auctions::SecondPrice(std::vector<numericaldists::Distribution>{
      boost::math::uniform_distribution<>(20, 40),
      boost::math::uniform_distribution<>(50, 70)});
  epsilon = 0.015;
  auto bid_func = [](float x) { return x; };
  auto bid_func2 = [](float x) { return x; };
  auction.AcceptStrategy(bid_func, 0);
  auction.AcceptStrategy(bid_func2, 1);
  float fit = auction.GetFitness(bid_func, 0);
  EXPECT_NEAR(0, fit, epsilon);
  float fit2 = auction.GetFitness(bid_func2, 1);
  EXPECT_NEAR(30, fit2, epsilon);
}

}  // namespace gatests
