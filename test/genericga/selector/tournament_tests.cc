#include "genericga/selector/tournament.h"

#include "genericga/vector_ops.h"

#include <gtest/gtest.h>

#include <memory>
#include <vector>

namespace gatests {

class TournamentSelectorTest : public ::testing::Test {
 public:
  TournamentSelectorTest() {}

 protected:
  virtual void SetUp() {
    sel = std::make_unique<genericga::selector::Tournament>(3);
  }
  std::vector<float> fitnesses = {1, 3, 6, -2};
  std::vector<int> counts = {1, 1, 1, 1};
  std::unique_ptr<genericga::selector::Tournament> sel;
};

TEST_F(TournamentSelectorTest, TournamentRoundTest) {
  auto ranks = genericga::GetRankingsWithTies(fitnesses, genericga::MinRank);
  ASSERT_EQ(2, sel->TournamentRound(ranks, std::vector<int>{0, 2, 2}));
  ASSERT_EQ(1, sel->TournamentRound(ranks, std::vector<int>{0, 3, 1, 0}));
}

}  // namespace gatests
