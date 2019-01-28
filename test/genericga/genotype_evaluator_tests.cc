#include <gtest/gtest.h>

#include <vector>

#include "genericga/genotype_evaluator.h"
#include "genericga/phenotype_strategy.h"
#include "sample_fitness_calculator.h"
#include "sample_phenotype_converter.h"

namespace gatests {

class GenotypeEvaluatorTest : public ::testing::Test {
 public:
  GenotypeEvaluatorTest()
      : gene_eval(SamplePhenotypeConverter(100), SampleFitnessCalculator(5)) {}

 protected:
  virtual void SetUp() {}
  genericga::GenotypeEvaluator<int, int> gene_eval;
};

TEST_F(GenotypeEvaluatorTest, GetFitnessTest) {
  EXPECT_EQ(24, gene_eval.GetFitness(9));
}

TEST_F(GenotypeEvaluatorTest, GetFitnessesTest) {
  std::vector<float> exp_fits{24, 80, 69};
  EXPECT_EQ(exp_fits, gene_eval.GetFitnesses(std::vector<int>{9, 5, 6}));
}

TEST_F(GenotypeEvaluatorTest, GetPhenotypeTest) {
  EXPECT_EQ(19, gene_eval.GetPhenotype(9));
}

TEST_F(GenotypeEvaluatorTest, GetPhenotypeFitnessTest) {
  EXPECT_EQ(14, gene_eval.GetPhenotypeFitness(9));
}

TEST_F(GenotypeEvaluatorTest, GetPhenotypeStrategyTest) {
  genericga::PhenotypeStrategy<int> strat{19, 24};
  auto res = gene_eval.GetPhenotypeStrategy(9);
  EXPECT_EQ(strat.phenotype, res.phenotype);
  EXPECT_EQ(strat.fitness, res.fitness);
}

}  // namespace gatests
