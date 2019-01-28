#include "genericga/selector/tournament.h"
#include "../sample_fitness_collection.h"

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
  SampleFitnessCollection pop;
  std::unique_ptr<genericga::selector::Tournament> sel;
};

TEST_F(TournamentSelectorTest, TournamentRoundTest) {
  ASSERT_EQ(2, sel->TournamentRound(pop.GetFitnessRankings(),
                                    std::vector<int>{0, 2, 2}));
  ASSERT_EQ(1, sel->TournamentRound(pop.GetFitnessRankings(),
                                    std::vector<int>{0, 3, 1, 0}));
}

}  // namespace gatests
