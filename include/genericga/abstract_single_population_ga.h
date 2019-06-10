#ifndef GENERICGA_ABSTRACT_SINGLE_POPULATION_GA_H_
#define GENERICGA_ABSTRACT_SINGLE_POPULATION_GA_H_

#include <omp.h>
#include <memory>
#include <vector>

#include "genericga/children_factory.h"
#include "genericga/phenotype_strategy.h"
#include "genericga/population.h"
#include "genericga/selector.h"

namespace genericga {

template <class Phen>
class AbstractSinglePopulationGA {
 public:
  virtual void RunRound(int n) = 0;
  virtual void SetFitnessCalculator(
      std::function<std::vector<float>(const std::vector<Phen>&)> fit_calc) = 0;
  virtual void SetFitnessCalculator(
      std::function<float(const Phen&)> fit_calc) {
    auto par_fit = [fit_calc](const std::vector<Phen>& phens) {
      std::vector<float> fits(phens.size());
#pragma omp parallel for
      for (int i = 0; i < phens.size(); ++i) {
        fits[i] = fit_calc(phens[i]);
      }
      return fits;
    };
    SetFitnessCalculator(par_fit);
  }
  virtual PhenotypeStrategy<Phen> SelectStrategy(Selector& sel) = 0;
  virtual std::vector<PhenotypeStrategy<Phen>> SelectStrategies(Selector& sel,
                                                                int n) = 0;
  virtual PhenotypeStrategy<Phen> GetBestStrategy() = 0;
  virtual std::vector<PhenotypeStrategy<Phen>> GetBestStrategies(int n) = 0;
  virtual PhenotypeStrategy<Phen> GetCommonestStrategy() = 0;
  virtual std::vector<PhenotypeStrategy<Phen>> GetCommonestStrategies(int n) = 0;
  virtual ~AbstractSinglePopulationGA() {}
};

}  // namespace genericga

#endif  // GENERICGA_ABSTRACT_SINGLE_POPULATION_GA_H_
