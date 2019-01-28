#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "sample_fitness_collection.h"
#include "genericga/selector/keep_best.h"

namespace gatests {

class FitnessCollectionTest : public ::testing::Test {
 public:
  FitnessCollectionTest() {}

 protected:
  virtual void SetUp() {}

  SampleFitnessCollection fit;
  genericga::selector::KeepBest sel;
};

TEST_F(FitnessCollectionTest, OrderingsTest) {
  std::vector<int> orderings = {3, 0, 1, 2};
  ASSERT_EQ(orderings, fit.GetFitnessOrderings());
}

TEST_F(FitnessCollectionTest, RankingsTest) {
  std::vector<float> rankings = {1, 2, 3, 0};
  ASSERT_EQ(rankings, fit.GetFitnessRankings());
}

TEST_F(FitnessCollectionTest, SizeTest) {
  ASSERT_EQ(4, fit.Size());
}

}  // namespace gatests
