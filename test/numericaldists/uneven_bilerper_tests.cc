#include <gtest/gtest.h>

#include <vector>

#include "numericaldists/interval.h"
#include "numericaldists/uneven_bilerper.h"

namespace gatests {

using namespace numericaldists;

class UnevenBilerperTest : public ::testing::Test {
 public:
  UnevenBilerperTest() {}

 protected:
  virtual void SetUp() {}
  UnevenBilerper func =
      UnevenBilerper(std::vector<float>{5, 10}, std::vector<float>{0, 10},
                     std::vector<std::vector<float>>{std::vector<float>{4, 5},
                                                     std::vector<float>{1, 3}});
};

TEST_F(UnevenBilerperTest, GetBidInterior) {
  EXPECT_FLOAT_EQ(3.88, func(7, 8));
}

TEST_F(UnevenBilerperTest, GetBidExterior) {
  EXPECT_FLOAT_EQ(1, func(4, -1));
  EXPECT_FLOAT_EQ(2, func(7.5, -1));
  EXPECT_FLOAT_EQ(3, func(12, -1));
  EXPECT_FLOAT_EQ(2.5, func(4, 5));
  EXPECT_FLOAT_EQ(4, func(12, 5));
  EXPECT_FLOAT_EQ(4, func(4, 12));
  EXPECT_FLOAT_EQ(4.5, func(7.5, 12));
  EXPECT_FLOAT_EQ(5, func(12, 12));
}

TEST_F(UnevenBilerperTest, GetBidBoundary) {
  EXPECT_FLOAT_EQ(1, func(5, 0));
  EXPECT_FLOAT_EQ(4, func(5, 10));
  EXPECT_FLOAT_EQ(3, func(10, 0));
  EXPECT_FLOAT_EQ(5, func(10, 10));
}

}  // namespace gatests
