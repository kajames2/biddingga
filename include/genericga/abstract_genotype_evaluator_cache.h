#ifndef _GENERICGA_ABSTRACT_GENOTYPE_EVALUATOR_CACHE_H_
#define _GENERICGA_ABSTRACT_GENOTYPE_EVALUATOR_CACHE_H_

#include <vector>

#include "genericga/abstract_genotype_evaluator.h"
#include "genericga/phenotype_strategy.h"

namespace genericga {

template <class Gen, class Phen>
class AbstractGenotypeEvaluatorCache : public AbstractGenotypeEvaluator<Gen> {
 public:
  // Calculates the fitness value for genotypes passed in and stores them.
  virtual void SetCache(const std::vector<Gen>& genotypes) = 0;
  virtual void SetFitnessCalculator(
      std::function<float(const Phen&)> fit_calc) = 0;
  virtual PhenotypeStrategy<Phen> GetPhenotypeStrategy(
      const Gen& gene) const = 0;
  virtual std::vector<PhenotypeStrategy<Phen>> GetPhenotypeStrategies(
      const std::vector<Gen>& genes) const = 0;
};

}  // namespace genericga
#endif  // _GENERICGA_ABSTRACT_GENOTYPE_EVALUATOR_CACHE_H_3
