#ifndef GENERICGA_COMPOSITE_GA_H_
#define GENERICGA_COMPOSITE_GA_H_

#include <memory>
#include <vector>

#include "genericga/abstract_single_population_ga.h"
#include "genericga/phenotype_strategy.h"
#include "genericga/selector.h"
#include "genericga/selector/keep_best.h"
#include "genericga/selector/keep_commonest.h"

namespace genericga {

template <class Phen>
class CompositeGA : public AbstractSinglePopulationGA<Phen> {
 public:
  CompositeGA(
      std::vector<std::shared_ptr<AbstractSinglePopulationGA<Phen>>> gas,
      std::function<
          PhenotypeStrategy<Phen>(const std::vector<PhenotypeStrategy<Phen>>&)>
          strat_combiner)
      : gas_(gas),
        strat_combiner_(strat_combiner),
        best_sel_(),
        commonest_sel_() {}

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

    return comb_strats;
  }

  PhenotypeStrategy<Phen> GetBestStrategy() override {
    return SelectStrategy(best_sel_);
  }

  std::vector<PhenotypeStrategy<Phen>> GetBestStrategies(int n) override {
    return SelectStrategies(best_sel_, n);
  }

  PhenotypeStrategy<Phen> GetCommonestStrategy() override {
    return SelectStrategy(commonest_sel_);
  }

  std::vector<PhenotypeStrategy<Phen>> GetCommonestStrategies(int n) override {
    return SelectStrategies(commonest_sel_, n);
  }

 private:
  std::vector<std::shared_ptr<AbstractSinglePopulationGA<Phen>>> gas_;
  std::function<PhenotypeStrategy<Phen>(
      const std::vector<PhenotypeStrategy<Phen>>&)>
      strat_combiner_;
  selector::KeepBest best_sel_;
  selector::KeepCommonest commonest_sel_;
};

}  // namespace genericga

#endif  // GENERICGA_COMPOSITE_GA_H_
