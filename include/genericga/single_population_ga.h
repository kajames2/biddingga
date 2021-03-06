#ifndef GENERICGA_SINGLE_POPULATION_GA_H_
#define GENERICGA_SINGLE_POPULATION_GA_H_

#include <memory>
#include <utility>
#include <vector>

#include "genericga/abstract_single_population_ga.h"
#include "genericga/children_factory.h"
#include "genericga/phenotype_strategy.h"
#include "genericga/population.h"
#include "genericga/selector.h"
#include "genericga/selector/elitism_decorator.h"
#include "genericga/selector/keep_best.h"
#include "genericga/selector/keep_commonest.h"
#include "genericga/selector/ranked_weighted.h"

namespace genericga {

template <class Gen, class Phen>
class SinglePopulationGA : public AbstractSinglePopulationGA<Phen> {
 public:
  SinglePopulationGA(Population<Gen, Phen> init_pop,
                     std::unique_ptr<ChildrenFactory<Gen>> children_fact,
                     std::unique_ptr<Selector> survivor_selector =
                         std::make_unique<selector::ElitismDecorator>(
                             std::make_unique<selector::RankedWeighted>(0.6),
                             1));
  void RunRound(int n = 1) override;
  void SetFitnessCalculator(
      std::function<std::vector<float>(const std::vector<Phen>&)> fit_calc)
      override;

  std::vector<PhenotypeStrategy<Phen>> GetPopulation() const {
    return pop_.GetPhenotypeStrategies();
  }

  PhenotypeStrategy<Phen> SelectStrategy(Selector& sel) override {
    return pop_.SelectPhenotypeStrategy(sel);
  }

  std::vector<PhenotypeStrategy<Phen>> SelectStrategies(Selector& sel,
                                                        int n) override {
    return pop_.SelectPhenotypeStrategies(sel, n);
  }

  PhenotypeStrategy<Phen> GetBestStrategy() override {
    return SelectStrategy(best_sel);
  }

  std::vector<PhenotypeStrategy<Phen>> GetBestStrategies(int n) override {
    return SelectStrategies(best_sel, n);
  }

  PhenotypeStrategy<Phen> GetCommonestStrategy() override {
    return SelectStrategy(commonest_sel);
  }

  std::vector<PhenotypeStrategy<Phen>> GetCommonestStrategies(int n) override {
    return SelectStrategies(commonest_sel, n);
  }

 private:
  void RunSingleRound();

  Population<Gen, Phen> pop_;
  int n_strategies_;
  int n_children_;
  selector::KeepBest best_sel;
  selector::KeepCommonest commonest_sel;
  std::unique_ptr<Selector> survivor_selector_;
  std::unique_ptr<ChildrenFactory<Gen>> children_fact_;
};

template <class Gen, class Phen>
SinglePopulationGA<Gen, Phen>::SinglePopulationGA(
    Population<Gen, Phen> init_pop,
    std::unique_ptr<ChildrenFactory<Gen>> children_fact,
    std::unique_ptr<Selector> survivor_selector)
    : pop_(std::move(init_pop)),
      n_strategies_(pop_.GetTotalCount()),
      n_children_(pop_.GetTotalCount()),
      survivor_selector_(std::move(survivor_selector)),
      children_fact_(std::move(children_fact)),
      best_sel(),
      commonest_sel() {}

template <class Gen, class Phen>
void SinglePopulationGA<Gen, Phen>::RunSingleRound() {
  auto children = children_fact_->GetChildren(pop_, n_children_);
  pop_.AddGenotypes(std::move(children));
  pop_.Survival(*survivor_selector_, n_strategies_);
}

template <class Gen, class Phen>
void SinglePopulationGA<Gen, Phen>::RunRound(int n) {
  for (int i = 0; i < n; ++i) {
    RunSingleRound();
  }
}

template <class Gen, class Phen>
void SinglePopulationGA<Gen, Phen>::SetFitnessCalculator(
    std::function<std::vector<float>(const std::vector<Phen>&)> fit_calc) {
  pop_.SetFitnessCalculator(fit_calc);
}

}  // namespace genericga

#endif  // GENERICGA_SINGLE_POPULATION_GA_H_
