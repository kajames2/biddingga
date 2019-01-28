#include <gtest/gtest.h>

#include <vector>

#include "genericga/population.h"
#include "genericga/selector/keep_best.h"
#include "sample_fitness_calculator.h"
#include "sample_phenotype_converter.h"

namespace gatests {

class PopulationTest : public ::testing::Test {
 public:
  PopulationTest()
      : pop(SamplePhenotypeConverter(100), SampleFitnessCalculator(5),
            std::vector<int>{5, 1, 9}) {}

 protected:
  virtual void SetUp() {}
  genericga::Population<int, int> pop;
  genericga::selector::KeepBest sel;
};

typedef PopulationTest PopulationDeathTest;

TEST_F(PopulationTest, GetFitnessesTest) {
  std::vector<float> fit = {80, 104, 24};
  EXPECT_EQ(fit, pop.GetFitnesses());
}

TEST_F(PopulationTest, SetFitnessUpdates) {
  std::vector<float> fit1 = {80, 104, 24};
  std::vector<float> fit2 = {81, 105, 25};
  EXPECT_EQ(fit1, pop.GetFitnesses());
  pop.SetFitnessCalculator(SampleFitnessCalculator(6));
  EXPECT_EQ(fit2, pop.GetFitnesses());
}

TEST_F(PopulationTest, AddGenotypesTest) {
  std::vector<float> exp_fits{80, 104, 24, 101, 101};
  pop.AddGenotypes(std::vector<int>(2, 2));
  EXPECT_EQ(exp_fits, pop.GetFitnesses());
}

TEST_F(PopulationTest, SetGenotypesTest) {
  std::vector<float> exp_fits{101, 96};
  pop.SetGenotypes(std::vector<int>{2, 3});
  EXPECT_EQ(exp_fits, pop.GetFitnesses());
}

TEST_F(PopulationTest, GetAllPhenotypeStrategiesTest) {
  genericga::PhenotypeStrategy<int> strat1{75, 80};
  genericga::PhenotypeStrategy<int> strat2{99, 104};
  genericga::PhenotypeStrategy<int> strat3{19, 24};
  std::vector<genericga::PhenotypeStrategy<int>> phens = {strat1, strat2,
                                                          strat3};
  auto res = pop.GetAllPhenotypeStrategies();
  EXPECT_EQ(phens[0].phenotype, res[0].phenotype);
  EXPECT_EQ(phens[0].fitness, res[0].fitness);
  EXPECT_EQ(phens[2].phenotype, res[2].phenotype);
  EXPECT_EQ(phens[2].fitness, res[2].fitness);
  EXPECT_EQ(phens.size(), res.size());
}

TEST_F(PopulationTest, GetSelectPhenotypeStrategyTest) {
  genericga::PhenotypeStrategy<int> strat{99, 104};
  auto res = pop.SelectPhenotypeStrategy(sel);
  EXPECT_EQ(strat.phenotype, res.phenotype);
  EXPECT_EQ(strat.fitness, res.fitness);
}

TEST_F(PopulationDeathTest, SelectWhenEmptyDeathTest) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  pop.SetGenotypes(std::vector<int>());
  EXPECT_DEATH(pop.SelectPhenotypeStrategy(sel), "");
}

TEST_F(PopulationTest, SelectPhenotypeStrategiesTest) {
  genericga::PhenotypeStrategy<int> strat1{99, 104};
  genericga::PhenotypeStrategy<int> strat2{75, 80};
  std::vector<genericga::PhenotypeStrategy<int>> phens = {strat1, strat2};
  auto res = pop.SelectPhenotypeStrategies(sel, 2);
  EXPECT_EQ(phens[0].phenotype, res[0].phenotype);
  EXPECT_EQ(phens[0].fitness, res[0].fitness);
  EXPECT_EQ(phens[1].phenotype, res[1].phenotype);
  EXPECT_EQ(phens[1].fitness, res[1].fitness);
  EXPECT_EQ(phens.size(), res.size());
}

}  // namespace gatests
