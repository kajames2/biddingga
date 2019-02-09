#ifndef _GENERICGA_MULTIPOP_SUB_GA_ADAPTER_H_
#define _GENERICGA_MULTIPOP_SUB_GA_ADAPTER_H_

#include "genericga/abstract_single_population_ga.h"
#include "genericga/multipop/abstract_sub_ga.h"
#include "genericga/selector.h"
#include "genericga/selector/tournament.h"

#include <algorithm>
#include <functional>
#include <memory>

namespace genericga {
namespace multipop {

template <class Environment, class Phen>
class SubGAAdapter : public AbstractSubGA<Environment> {
 public:
  SubGAAdapter(std::unique_ptr<AbstractSinglePopulationGA<Phen>> ga, int id,
               int priority = 0,
               std::unique_ptr<Selector> selector =
                   std::make_unique<selector::Tournament>(3));
  void RunRound(const std::vector<Environment>& envs) override;
  void SubmitPlayStrat(Environment& env) override;
  PhenotypeStrategy<Phen> SelectStrategy(Selector& sel) {
    return ga_->SelectStrategy(sel);
  }
  std::vector<PhenotypeStrategy<Phen>> SelectStrategies(Selector& sel, int n) {
    return ga_->SelectStrategies(sel, n);
  }
  PhenotypeStrategy<Phen> GetBestStrategy() { return ga_->GetBestStrategy(); }
  std::vector<PhenotypeStrategy<Phen>> GetBestStrategies(int n) {
    return ga_->GetBestStrategies(n);
  }
  std::vector<PhenotypeStrategy<Phen>> GetPopulation() const {
    return ga_->GetPopulation();
  }

 private:
  std::function<std::vector<float>(const std::vector<Phen>&)>
  GetFitnessFunction(const std::vector<Environment>& envs) const;
  static float GetAccumulatedFitness(const std::vector<Environment>& envs,
                                     int id, const Phen& phen);

  std::unique_ptr<Selector> selector_;
  std::unique_ptr<AbstractSinglePopulationGA<Phen>> ga_;
};

template <class Environment, class Phen>
SubGAAdapter<Environment, Phen>::SubGAAdapter(
    std::unique_ptr<AbstractSinglePopulationGA<Phen>> ga, int id, int priority,
    std::unique_ptr<Selector> selector)
    : AbstractSubGA<Environment>(id, priority),
      selector_(std::move(selector)),
      ga_(std::move(ga)) {}

template <class Environment, class Phen>
void SubGAAdapter<Environment, Phen>::RunRound(
    const std::vector<Environment>& envs) {
  ga_->SetFitnessCalculator(GetFitnessFunction(envs));
  ga_->RunRound(1);
}

template <class Environment, class Phen>
void SubGAAdapter<Environment, Phen>::SubmitPlayStrat(Environment& env) {
  env.AcceptStrategy(ga_->SelectStrategy(*selector_).phenotype, this->GetID());
}

template <class Environment, class Phen>
std::function<std::vector<float>(const std::vector<Phen>&)>
SubGAAdapter<Environment, Phen>::GetFitnessFunction(
    const std::vector<Environment>& envs) const {
  int id = this->GetID();
  return [&envs, id](const std::vector<Phen>& phens) -> std::vector<float> {
    std::vector<float> tots(phens.size(), 0.0);
    for (const auto& env : envs) {
      auto env_fits = env.GetFitness(phens, id);
      std::transform(tots.begin(), tots.end(), env_fits.begin(), tots.begin(),
                     std::plus<float>());
    }
    for (auto& tot : tots) {
      tot /= envs.size();
    }
    return tots;
  };
}

}  // namespace multipop
}  // namespace genericga

#endif  // _GENERICGA_MULTIPOP_SINGLE_POPULATION_GA_TO_SUB_GA_ADAPTER_H_
