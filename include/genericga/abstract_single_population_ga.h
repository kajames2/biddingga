#ifndef _GENERICGA_ABSTRACT_SINGLE_POPULATION_GA_H_
#define _GENERICGA_ABSTRACT_SINGLE_POPULATION_GA_H_

#include <memory>
#include <vector>

#include "genericga/children_factory.h"
#include "genericga/phenotype_strategy.h"
#include "genericga/population.h"
#include "genericga/selector.h"
#include "genericga/selector/ranked_weighted.h"

namespace genericga {

template <class Phen>
class AbstractSinglePopulationGA {
 public:
  virtual void RunRound(int n) = 0;
  virtual void SetFitnessCalculator(std::function<float(const Phen&)> fit_calc) = 0;
  virtual PhenotypeStrategy<Phen> SelectStrategy(Selector& sel) = 0;
  virtual std::vector<PhenotypeStrategy<Phen>> SelectStrategies(Selector& sel, int n) = 0;
  virtual PhenotypeStrategy<Phen> GetBestStrategy() = 0;
  virtual std::vector<PhenotypeStrategy<Phen>> GetBestStrategies(int n) = 0;
  virtual std::vector<PhenotypeStrategy<Phen>> GetPopulation() const = 0;
  virtual ~AbstractSinglePopulationGA() {}
};

}  // namespace genericga

#endif  // _GENERICGA_ABSTRACT_SINGLE_POPULATION_GA_H_
