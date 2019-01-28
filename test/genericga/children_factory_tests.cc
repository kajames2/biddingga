#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "genericga/children_factory.h"
#include "genericga/selector/keep_best.h"
#include "sample_crossover.h"
#include "sample_mutation.h"
#include "sample_genotype_population.h"

namespace gatests {

class ChildrenFactoryTest : public ::testing::Test {
 public:
  ChildrenFactoryTest() {}

 protected:
  virtual void SetUp() {
    auto cross = std::make_unique<SampleCrossover>();
    auto mut = std::make_unique<SampleMutation>();
    auto sel = std::make_unique<genericga::selector::KeepBest>();
    fact = std::make_unique<genericga::ChildrenFactory<int>>(
        std::move(cross), std::move(mut), std::move(sel));
  }

  std::unique_ptr<genericga::ChildrenFactory<int>> fact;
  SampleGenotypePopulation pop;
};

TEST_F(ChildrenFactoryTest, GetChildrenTest) {
  auto children = fact->GetChildren(pop, 2);
  ASSERT_EQ(2, children.size());
  ASSERT_EQ(48, children[0]);
  ASSERT_EQ(48, children[1]);
}

}  // namespace gatests
