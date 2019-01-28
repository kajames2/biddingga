#include <gtest/gtest.h>

#include <vector>

#include "genericga/genotype_evaluator_single_cache.h"
#include "genericga/phenotype_strategy.h"
#include "sample_fitness_calculator.h"
#include "sample_phenotype_converter.h"

namespace gatests {

class GenotypeEvaluatorSingleCacheTest : public ::testing::Test {
 public:
  GenotypeEvaluatorSingleCacheTest()
      : gene_eval(SamplePhenotypeConverter(100), SampleFitnessCalculator(5)) {}

 protected:
  virtual void SetUp() {}
  genericga::GenotypeEvaluatorSingleCache<int, int> gene_eval;
};

TEST_F(GenotypeEvaluatorSingleCacheTest, GetFitnessTest) {
  EXPECT_EQ(24, gene_eval.GetFitness(9));
}

TEST_F(GenotypeEvaluatorSingleCacheTest, SetFitnessUpdatesCacheTest) {
  gene_eval.GetFitness(9);
  gene_eval.SetFitnessCalculator(SampleFitnessCalculator(6));
  EXPECT_EQ(25, gene_eval.GetFitness(9));
}

TEST_F(GenotypeEvaluatorSingleCacheTest, GetFitnessesTest) {
  std::vector<float> exp_fits{24, 80, 69};
  EXPECT_EQ(exp_fits, gene_eval.GetFitnesses(std::vector<int>{9, 5, 6}));
}

TEST_F(GenotypeEvaluatorSingleCacheTest, GetPhenotypeStrategyTest) {
  genericga::PhenotypeStrategy<int> strat{19, 24};
  auto res = gene_eval.GetPhenotypeStrategy(9);
  EXPECT_EQ(strat.phenotype, res.phenotype);
  EXPECT_EQ(strat.fitness, res.fitness);
}

TEST_F(GenotypeEvaluatorSingleCacheTest, GetPhenotypeStrategiesTest) {
  genericga::PhenotypeStrategy<int> strat1{19, 24};
  genericga::PhenotypeStrategy<int> strat2{0, 5};
  std::vector<int> genes = {9,10};
  std::vector<genericga::PhenotypeStrategy<int>> phens = {strat1, strat2};
  auto res = gene_eval.GetPhenotypeStrategies(genes);
  EXPECT_EQ(phens[0].phenotype, res[0].phenotype);
  EXPECT_EQ(phens[0].fitness, res[0].fitness);
  EXPECT_EQ(phens[1].phenotype, res[1].phenotype);
  EXPECT_EQ(phens[1].fitness, res[1].fitness);
  EXPECT_EQ(phens.size(), res.size());
}

}  // namespace gatests
