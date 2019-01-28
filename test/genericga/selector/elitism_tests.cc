#include "genericga/selector/elitism_decorator.h"
#include "genericga/selector/keep_best.h"
#include "../sample_fitness_collection.h"

#include <gtest/gtest.h>
#include <map>
#include <memory>

namespace gatests {

class ElitismTest : public ::testing::Test {
 public:
  ElitismTest() {}

 protected:
  virtual void SetUp() {
    auto inner_sel = std::make_unique<genericga::selector::KeepBest>();
    sel = std::make_unique<genericga::selector::ElitismDecorator>(
        std::move(inner_sel), 2);
  }

  std::unique_ptr<genericga::selector::ElitismDecorator> sel;
  SampleFitnessCollection pop;
};

TEST_F(ElitismTest, SelectIndicesTest) {
  auto vec = sel->SelectIndices(pop, 5);
  ASSERT_EQ(5, vec.size());
  std::map<int, int> counts;
  for (int val : vec) {
    if (counts.find(val) == counts.end()) {
      counts.emplace(val, 1);
    } else {
      counts[val] += 1;
    }
  }
  ASSERT_EQ(1, counts[0]);
  ASSERT_EQ(2, counts[1]);
  ASSERT_EQ(2, counts[2]);
}

TEST_F(ElitismTest, AllElitesTest) {
  auto vec = sel->SelectIndices(pop, 1);
  ASSERT_EQ(1, vec.size());
  ASSERT_EQ(2, vec[0]);
}

}  // namespace gatests
