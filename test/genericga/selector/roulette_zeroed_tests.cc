#include "genericga/selector/roulette_zeroed.h"

#include <gtest/gtest.h>

#include <vector>
#include <memory>

namespace gatests {

class RouletteZeroedTest : public ::testing::Test {
public:
  RouletteZeroedTest() {}

protected:
  virtual void SetUp() {
    sel = std::make_unique<genericga::selector::RouletteZeroed>();
  }
  std::vector<float> fitnesses = {1, 3, 6, -2};
  std::unique_ptr<genericga::selector::RouletteZeroed> sel;
};

TEST_F(RouletteZeroedTest, CalculateWeightsTest) {
  auto vec = sel->CalculateWeights(fitnesses);
  ASSERT_EQ(4, vec.size());
  ASSERT_DOUBLE_EQ(3.0, vec[0]);
  ASSERT_EQ(5.0, vec[1]);
  ASSERT_EQ(0.0, vec[3]);
}

} // namespace gatests
