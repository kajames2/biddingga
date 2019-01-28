#include "genericga/selector/roulette_simple.h"
#include "../sample_fitness_collection.h"

#include <gtest/gtest.h>

#include <memory>
#include <vector>

namespace gatests {

class RouletteSimpleTest : public ::testing::Test {
 public:
  RouletteSimpleTest() {}

 protected:
  virtual void SetUp() {
    sel = std::make_unique<genericga::selector::RouletteSimple>();
  }

  SampleFitnessCollection pop;
  std::unique_ptr<genericga::selector::RouletteSimple> sel;
};

TEST_F(RouletteSimpleTest, CalculateWeightsTest) {
  auto vec = sel->CalculateWeights(pop);
  ASSERT_EQ(4, vec.size());
  ASSERT_DOUBLE_EQ(1.0, vec[0]);
  ASSERT_EQ(3.0, vec[1]);
  ASSERT_EQ(-2.0, vec[3]);
}

}  // namespace gatests
