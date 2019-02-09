#include "genericga/selector/ranked_weighted.h"

#include <gtest/gtest.h>

#include <vector>
#include <memory>

namespace gatests {

class RankedWeightedTest : public ::testing::Test {
public:
  RankedWeightedTest() {}

protected:
  virtual void SetUp() {
    sel = std::make_unique<genericga::selector::RankedWeighted>(0.6);
  }
  std::vector<float> fitnesses = {1, 3, 6, -2};
  std::unique_ptr<genericga::selector::RankedWeighted> sel;
};

TEST_F(RankedWeightedTest, CalculateWeightsTest) {
  auto vec = sel->CalculateWeights(fitnesses);
  EXPECT_EQ(4, vec.size());
  EXPECT_FLOAT_EQ(0.2, vec[0]);
  EXPECT_FLOAT_EQ(0.3, vec[1]);
  EXPECT_FLOAT_EQ(0.4, vec[2]);
  EXPECT_FLOAT_EQ(0.1, vec[3]);
}

} // namespace gatests
