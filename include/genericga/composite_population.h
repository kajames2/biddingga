#ifndef GENERICGA_COMPOSITE_GA_H_
#define GENERICGA_COMPOSITE_GA_H_

#include <memory>
#include <vector>
#include <omp.h>

#include "genericga/abstract_single_population_ga.h"
#include "genericga/phenotype_strategy.h"
#include "genericga/selector.h"
#include "genericga/selector/keep_best.h"

namespace genericga {

template <class Phen>
class CompositeGA : public AbstractSinglePopulationGA<Phen> {
 public:
  void RunRound(int n) override {
    for (auto& ga : gas_) {
      ga->RunRound(n);
    }
  };
  void SetFitnessCalculator(
      std::function<std::vector<float>(const std::vector<Phen>&)> fit_calc)
      override {
    for (auto& ga : gas_) {
      ga->SetFitnessCalculator(fit_calc);
    }
  }

  PhenotypeStrategy<Phen> SelectStrategy(Selector& sel) override {
    std::vector<PhenotypeStrategy<Phen>> strats;
    for (auto& ga : gas_) {
      strats.push_back(ga->SelectStrategy(sel));
    }
    return strat_combiner_(strats);
  };

  std::vector<PhenotypeStrategy<Phen>> SelectStrategies(Selector& sel,
                                                        int n) override {
    std::vector<std::vector<PhenotypeStrategy<Phen>>> strat_sets(n);
    for (auto& ga : gas_) {
      auto strats = ga->SelectStrategies(sel, n);
      for (int i = 0; i < n; ++i) {
        strat_sets[i].push_back(strats[i]);
      }
    }
    std::vector<PhenotypeStrategy<Phen>> comb_strats;
    std::transform(strat_sets.begin(), strat_sets.end(),
                   std::back_inserter(comb_strats),
                   [this](const std::vector<PhenotypeStrategy<Phen>>& strats) {
                     return strat_combiner_(strats);
                   });
  }

  PhenotypeStrategy<Phen> GetBestStrategy() override {
    return SelectStrategy(selector::KeepBest());
  }

  std::vector<PhenotypeStrategy<Phen>> GetBestStrategies(int n) override {
    return SelectStrategies(selector::KeepBest(), n);
  }

 private:
  std::vector<std::shared_ptr<AbstractSinglePopulationGA<Phen>>> gas_;
  std::function<PhenotypeStrategy<Phen>(
      const std::vector<PhenotypeStrategy<Phen>>&)>
      strat_combiner_;
};

}  // namespace genericga

#endif  // GENERICGA_COMPOSITE_GA_H_
