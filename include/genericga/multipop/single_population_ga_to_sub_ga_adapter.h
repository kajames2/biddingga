#ifndef _GENERICGA_MULTIPOP_SINGLE_POPULATION_GA_TO_SUB_GA_ADAPTER_H_
#define _GENERICGA_MULTIPOP_SINGLE_POPULATION_GA_TO_SUB_GA_ADAPTER_H_

#include "genericga/abstract_single_population_ga.h"
#include "genericga/multipop/abstract_sub_ga.h"
#include "genericga/selector.h"
#include "genericga/selector/tournament.h"

#include <memory>

namespace genericga {
namespace multipop {

template <class Environment, class Phen>
class SinglePopulationGAToSubGAAdapter : public AbstractSubGA<Environment> {
 public:
  SinglePopulationGAToSubGAAdapter(
      std::unique_ptr<AbstractSinglePopulationGA<Phen>> ga, int id,
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
  std::function<float(const Phen&)> GetFitnessFunction(
      const std::vector<Environment>& envs) const;
  static float GetAccumulatedFitness(const std::vector<Environment>& envs,
                                     int id, const Phen& phen);

  std::unique_ptr<Selector> selector_;
  std::unique_ptr<AbstractSinglePopulationGA<Phen>> ga_;
};

template <class Environment, class Phen>
SinglePopulationGAToSubGAAdapter<Environment, Phen>::
    SinglePopulationGAToSubGAAdapter(
        std::unique_ptr<AbstractSinglePopulationGA<Phen>> ga, int id,
        int priority, std::unique_ptr<Selector> selector)
    : AbstractSubGA<Environment>(id, priority),
      selector_(std::move(selector)),
      ga_(std::move(ga)) {}

template <class Environment, class Phen>
void SinglePopulationGAToSubGAAdapter<Environment, Phen>::RunRound(
    const std::vector<Environment>& envs) {
  ga_->SetFitnessCalculator(GetFitnessFunction(envs));
  ga_->RunRound(1);
}

template <class Environment, class Phen>
void SinglePopulationGAToSubGAAdapter<Environment, Phen>::SubmitPlayStrat(
    Environment& env) {
  env.AcceptStrategy(ga_->SelectStrategy(*selector_).phenotype, this->GetID());
}

template <class Environment, class Phen>
std::function<float(const Phen&)>
SinglePopulationGAToSubGAAdapter<Environment, Phen>::GetFitnessFunction(
    const std::vector<Environment>& envs) const {
  int id = this->GetID();
  return [&envs, id](const Phen& phen) -> float {
    float tot = 0.0;
    for (const auto& env : envs) {
      tot += env.GetFitness(phen,id);
    }
    return tot / envs.size();
  };
}

}  // namespace multipop
}  // namespace genericga

#endif  // _GENERICGA_MULTIPOP_SINGLE_POPULATION_GA_TO_SUB_GA_ADAPTER_H_
