#ifndef GENERICGA_MULTIPOP_SUB_GA_ADAPTER_H_
#define GENERICGA_MULTIPOP_SUB_GA_ADAPTER_H_

#include <algorithm>
#include <functional>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include "genericga/abstract_single_population_ga.h"
#include "genericga/multipop/abstract_sub_ga.h"
#include "genericga/selector.h"
#include "genericga/selector/keep_best.h"
#include "genericga/selector/tournament.h"

namespace genericga {
namespace multipop {

template <typename, typename T>
struct HasGetFitness {
  static_assert(std::integral_constant<T, false>::value,
                "Second template parameter needs to be of function type.");
};

template <typename C, typename Ret, typename... Args>
struct HasGetFitness<C, Ret(Args...)> {
 private:
  template <typename T>
  static constexpr auto check(T*) -> typename std::is_same<
      decltype(std::declval<T>().GetFitness(std::declval<Args>()...)),
      Ret>::type;

  template <typename>
  static constexpr std::false_type check(...);

  typedef decltype(check<C>(0)) type;

 public:
  static constexpr bool value = type::value;
};

template <class Environment, class Phen>
class SubGAAdapter : public AbstractSubGA<Environment> {
 public:
  SubGAAdapter(std::unique_ptr<AbstractSinglePopulationGA<Phen>> ga, int id,
               int priority = 0,
               std::unique_ptr<Selector> selector =
                   std::make_unique<selector::Tournament>(3));
  void RunRound(const std::vector<Environment*>& envs) override;
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
  PhenotypeStrategy<Phen> GetCommonestStrategy() {
    return ga_->GetCommonestStrategy();
  }
  std::vector<PhenotypeStrategy<Phen>> GetCommonestStrategies(int n) {
    return ga_->GetCommonestStrategies(n);
  }

 private:
  template <typename EnvMap>
  void SetFitnessCalculator(const std::vector<EnvMap*>& envs, std::true_type) {
    int id = this->GetID();
    auto func = [&envs,
                 id](const std::vector<Phen>& phens) -> std::vector<float> {
      std::vector<float> tots(phens.size(), 0.0);
      for (auto env : envs) {
        auto env_fits = env->GetFitness(phens, id);
        std::transform(tots.begin(), tots.end(), env_fits.begin(), tots.begin(),
                       std::plus<float>());
      }
      for (auto& tot : tots) {
        tot /= envs.size();
      }
      return tots;
    };
    ga_->SetFitnessCalculator(func);
  }

  template <typename EnvMap>
  void SetFitnessCalculator(const std::vector<EnvMap*>& envs, std::false_type) {
    int id = this->GetID();
    auto func = [&envs, id](const Phen& phen) -> float {
      float tot = 0.0;
      for (auto env : envs) {
        tot += env->GetFitness(phen, id);
      }
      return tot / envs.size();
    };
    ga_->SetFitnessCalculator(func);
  }

  template <typename EnvMap>
  void SetFitnessCalculator(const std::vector<EnvMap*>& envs) {
    SetFitnessCalculator(
        envs,
        std::integral_constant<
            bool,
            HasGetFitness<EnvMap, std::vector<float>(const std::vector<Phen>&,
                                                     int)>::value>());
  }

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
    const std::vector<Environment*>& envs) {
  SetFitnessCalculator(envs);
  ga_->RunRound(1);
}

template <class Environment, class Phen>
void SubGAAdapter<Environment, Phen>::SubmitPlayStrat(Environment& env) {
  env.AcceptStrategy(ga_->SelectStrategy(*selector_).phenotype, this->GetID());
}

}  // namespace multipop
}  // namespace genericga

#endif  // GENERICGA_MULTIPOP_SUB_GA_ADAPTER_H_
