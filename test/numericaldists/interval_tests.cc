#include <gtest/gtest.h>

#include <vector>

#include "numericaldists/interval.h"

namespace gatests {

using namespace numericaldists;

class IntervalTest : public ::testing::Test {
 public:
  IntervalTest() {}

 protected:
  virtual void SetUp() {}
};

TEST_F(IntervalTest, InIntervalTest) {
  auto interval = Interval{5, 10};
  EXPECT_FALSE(InInterval(interval, 0));
  EXPECT_TRUE(InInterval(interval, 5));
  EXPECT_TRUE(InInterval(interval, 8));
  EXPECT_FALSE(InInterval(interval, 10));
  EXPECT_FALSE(InInterval(interval, 15));
}

}  // namespace gatests
