#include <gtest/gtest.h>

#include "auctions/first_price.h"
#include "boost/math/distributions/uniform.hpp"
#include "numericaldists/distribution.h"
#include "numericaldists/function_ops.h"

#include <vector>

#include <eigen3/Eigen/Core>

using namespace Eigen;

namespace gatests {

class FirstPriceTest : public ::testing::Test {
 public:
  FirstPriceTest() {}

 protected:
  virtual void SetUp() {}

  ArrayXd values = ArrayXd::LinSpaced(101, 0, 1);
  auctions::FirstPrice auction = auctions::FirstPrice(
      std::vector<numericaldists::Distribution>{
          boost::math::uniform_distribution<>(0, 1),
          boost::math::uniform_distribution<>(0, 1)},
      101);
  float epsilon = 0.0001;
  numericaldists::Scatter bid_0 = {values, ArrayXd::Zero(1, 101)};
  numericaldists::Scatter bid_tiny = {values, ArrayXd::Ones(1, 101) * 0.01};
  numericaldists::Scatter bid_1 = {values, ArrayXd::Ones(1, 101)};
  numericaldists::Scatter bid_20 = {ArrayXd::LinSpaced(101, 20, 40), ArrayXd::Ones(1, 101) * 20};
  numericaldists::Scatter bid_202 = {ArrayXd::LinSpaced(101, 40, 60), ArrayXd::Ones(1, 101) * 20.2};
  numericaldists::Scatter bid_opt = {values, ArrayXd::LinSpaced(101, 0, 1) / 2};
  numericaldists::Scatter bid_neg = {values, ArrayXd::LinSpaced(101, 1, 0) / 2};
};

TEST_F(FirstPriceTest, AlwaysWinTest) {
  auction.AcceptStrategy(bid_0, 0);
  float fit = auction.GetFitness(bid_tiny, 1);
  EXPECT_NEAR(0.49, fit, epsilon);
}

// TEST_F(FirstPriceTest, AlwaysTieTest) {
//   auction.AcceptStrategy(bid_0, 0);
//   float fit = auction.GetFitness(bid_0, 1);
//   EXPECT_NEAR(0.25, fit, epsilon);
// }

TEST_F(FirstPriceTest, WinHalfTest) {
  auction.AcceptStrategy(bid_opt, 0);
  float fit = auction.GetFitness(bid_opt, 1);
  EXPECT_NEAR(.166666, fit, epsilon);
}

TEST_F(FirstPriceTest, DecreasingBidsTest) {
  auction.AcceptStrategy(bid_neg, 0);
  float fit = auction.GetFitness(bid_opt, 1);
  EXPECT_NEAR(.166666, fit, epsilon);
}

TEST_F(FirstPriceTest, DecreasingBids2Test) {
  auction.AcceptStrategy(bid_opt, 0);
  float fit = auction.GetFitness(bid_neg, 1);
  EXPECT_NEAR(0, fit, epsilon);
}

TEST_F(FirstPriceTest, NoIntersectDistsTest) {
  auction = auctions::FirstPrice(
      std::vector<numericaldists::Distribution>{
          boost::math::uniform_distribution<>(20, 40),
          boost::math::uniform_distribution<>(40, 60)},
      101);
   float epsilon = 0.001;
  auction.AcceptStrategy(bid_20, 0);
  auction.AcceptStrategy(bid_202, 1);
  float fit = auction.GetFitness(bid_20, 0);
  EXPECT_NEAR(0, fit, epsilon);
  float fit2 = auction.GetFitness(bid_202, 1);
  EXPECT_NEAR(29.8, fit2, epsilon);
}

TEST_F(FirstPriceTest, BidSegmentTest) {
  auction.AcceptStrategy(bid_opt, 0);
  numericaldists::Scatter bid_seg = {ArrayXd::LinSpaced(2, 0, 0.5), ArrayXd::LinSpaced(2, 1, 0.5)};
  float fit = auction.GetFitness(bid_seg, 1);
  EXPECT_NEAR(-0.25, fit, epsilon);
}

}  // namespace gatests
