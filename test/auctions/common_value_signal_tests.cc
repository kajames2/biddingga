#include <gtest/gtest.h>
#include <vector>
#include "auctions/common_value_signal.h"
#include "boost/math/distributions/uniform.hpp"
#include "numericaldists/distribution.h"
#include "numericaldists/function_ops.h"

namespace gatests {

class CommonValueSignalTest : public ::testing::Test {
 public:
  CommonValueSignalTest() {}

 protected:
  virtual void SetUp() {}
  auctions::CommonValueSignal auction = auctions::CommonValueSignal({2, 2}, 1, {-4, 3});
  float epsilon = 0.0001;
};

TEST_F(CommonValueSignalTest, AlwaysWinTest) {
  auction.AcceptStrategy([](float x) { return 3; }, 0);
  auction.AcceptStrategy([](float x) { return 1; }, 1);
  float fit = auction.GetFitness([](float x) { return 1; }, 1);
  EXPECT_NEAR(1, fit, epsilon);
}

TEST_F(CommonValueSignalTest, WinHalfTest) {
  auto bid_func = [](float x) { return 0; };
  auction.AcceptStrategy(bid_func, 0);
  auction.AcceptStrategy(bid_func, 1);
  float fit = auction.GetFitness(bid_func, 1);
  EXPECT_NEAR(-0.116666, fit, epsilon);
}

// TEST_F(CommonValueSignalTest, DecreasingBidsTest) {
//   auto bid_func = [](float x) { return 0.5 - x / 2; };
//   auction.AcceptStrategy(bid_func, 0);
//   float fit = auction.GetFitness([](float x) { return x / 2; }, 1);
//   EXPECT_NEAR(.166666, fit, epsilon);
// }

// TEST_F(CommonValueSignalTest, DecreasingBids2Test) {
//   auto bid_func = [](float x) { return x / 2; };
//   auction.AcceptStrategy(bid_func, 0);
//   float fit = auction.GetFitness([](float x) { return 0.5 - 0.5 * x; }, 1);
//   EXPECT_NEAR(0, fit, epsilon);
// }

// TEST_F(CommonValueSignalTest, NoIntersectDistsTest) {
//   auction =
//       auctions::CommonValueSignal(std::vector<numericaldists::Distribution>{
//           boost::math::uniform_distribution<>(20, 40),
//           boost::math::uniform_distribution<>(40, 60)});
//   float epsilon = 0.001;
//   auto bid_func = [](float x) { return 20; };
//   auction.AcceptStrategy(bid_func, 0);
//   auto bid_func2 = [](float x) { return 20.2; };
//   auction.AcceptStrategy(bid_func2, 1);
//   float fit = auction.GetFitness(bid_func, 0);
//   EXPECT_NEAR(0, fit, epsilon);
//   float fit2 = auction.GetFitness(bid_func2, 1);
//   EXPECT_NEAR(29.8, fit2, epsilon);
// }

}  // namespace gatests
