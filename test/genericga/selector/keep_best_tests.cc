#include "genericga/selector/keep_best.h"

#include <gtest/gtest.h>

#include <map>
#include <memory>

namespace gatests {

class KeepBestSelectorTest : public ::testing::Test {
 public:
  KeepBestSelectorTest() {}

 protected:
  virtual void SetUp() {
    sel = std::make_unique<genericga::selector::KeepBest>();
  }
  std::vector<float> fitnesses = {1, 3, 6, -2};
  std::vector<int> counts = {1, 1, 1, 1};
  std::unique_ptr<genericga::selector::KeepBest> sel;
};

TEST_F(KeepBestSelectorTest, SelectIndicesTest) {
  auto vec = sel->SelectIndices(fitnesses, counts, 2);
  ASSERT_EQ(2, vec.size());
  ASSERT_EQ(2, vec[0]);
  ASSERT_EQ(1, vec[1]);
}

}  // namespace gatests
