#include "genericga/selector/ranked_exponential.h"

#include <gtest/gtest.h>

#include <vector>
#include <memory>

namespace gatests {

class RankedExponentialTest : public ::testing::Test {
public:
  RankedExponentialTest() {}

protected:
  virtual void SetUp() {
    sel = std::make_unique<genericga::selector::RankedExponential>();
  }
  std::vector<float> fitnesses = {1, 3, 6, -2};
  std::unique_ptr<genericga::selector::RankedExponential> sel;
};

TEST_F(RankedExponentialTest, CalculateWeightsTest) {
  auto vec = sel->CalculateWeights(fitnesses);
  ASSERT_EQ(4, vec.size());
  ASSERT_FLOAT_EQ(1-std::exp(-1), vec[0]);
  ASSERT_FLOAT_EQ(1-std::exp(-3), vec[2]);
  ASSERT_FLOAT_EQ(1-std::exp(0), vec[3]);
}

} // namespace gatests
