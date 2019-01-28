#ifndef _GENERICGA_POPULATION_H_
#define _GENERICGA_POPULATION_H_

#include <algorithm>
#include <memory>
#include <vector>

#include "genericga/abstract_genotype_evaluator_cache.h"
#include "genericga/genotype_evaluator_double_cache.h"
#include "genericga/genotype_evaluator_single_cache.h"
#include "genericga/genotype_population.h"
#include "genericga/phenotype_strategy.h"

namespace genericga {

template <class Gen, class Phen>
class Population : public GenotypePopulation<Gen> {
 public:
  Population(
      std::unique_ptr<AbstractGenotypeEvaluatorCache<Gen, Phen>> gen_eval,
      std::vector<Gen> genes);

  Population(std::function<Phen(const Gen&)> phen_conv,
             std::function<float(const Phen&)> fit_calc,
             std::vector<Gen> genes);

  void SetFitnessCalculator(std::function<float(const Phen&)> fit_calc);
  std::vector<float> GetFitnesses() const override;

  void AddGenotypes(std::vector<Gen> genes) override;
  void SetGenotypes(std::vector<Gen> genes) override;

  PhenotypeStrategy<Phen> SelectPhenotypeStrategy(Selector& selector) const;
  std::vector<PhenotypeStrategy<Phen>> SelectPhenotypeStrategies(
      Selector& selector, int n) const;
  std::vector<PhenotypeStrategy<Phen>> GetAllPhenotypeStrategies() const;

  std::vector<Gen> GetAllGenotypes() const override { return genes_; }
  int Size() const override { return genes_.size(); }

 protected:
  Gen GetGenotype(int i) const override { return genes_[i]; }
  float GetFitness(int i) const override {
    return gen_eval_->GetFitness(genes_[i]);
  }
  PhenotypeStrategy<Phen> GetPhenotypeStrategy(int i) const {
    return gen_eval_->GetPhenotypeStrategy(GetGenotype(i));
  }

  std::vector<Gen> genes_;
  std::unique_ptr<AbstractGenotypeEvaluatorCache<Gen, Phen>> gen_eval_;
};

template <class Gen, class Phen>
Population<Gen, Phen>::Population(
    std::unique_ptr<AbstractGenotypeEvaluatorCache<Gen, Phen>> gen_eval,
    std::vector<Gen> genes)
    : genes_(std::move(genes)), gen_eval_(std::move(gen_eval)) {}

template <class Gen, class Phen>
Population<Gen, Phen>::Population(std::function<Phen(const Gen&)> phen_conv,
                                  std::function<float(const Phen&)> fit_calc,
                                  std::vector<Gen> genes)
    : genes_(std::move(genes)),
      gen_eval_(std::make_unique<GenotypeEvaluatorSingleCache<Gen, Phen>>(
          phen_conv, fit_calc)) {}

template <class Gen, class Phen>
void Population<Gen, Phen>::SetFitnessCalculator(
    std::function<float(const Phen&)> fit_calc) {
  gen_eval_->SetFitnessCalculator(std::move(fit_calc));
}

template <class Gen, class Phen>
std::vector<float> Population<Gen, Phen>::GetFitnesses() const {
  return gen_eval_->GetFitnesses(genes_);
}

template <class Gen, class Phen>
void Population<Gen, Phen>::AddGenotypes(std::vector<Gen> genes) {
  genes_.insert(genes_.end(), genes.begin(), genes.end());
  gen_eval_->SetCache(genes_);
}

template <class Gen, class Phen>
void Population<Gen, Phen>::SetGenotypes(std::vector<Gen> genes) {
  genes_ = genes;
  gen_eval_->SetCache(genes_);
}

template <class Gen, class Phen>
PhenotypeStrategy<Phen> Population<Gen, Phen>::SelectPhenotypeStrategy(
    Selector& selector) const {
  return gen_eval_->GetPhenotypeStrategy(
      GenotypePopulation<Gen>::SelectGenotype(selector));
}

template <class Gen, class Phen>
std::vector<PhenotypeStrategy<Phen>>
Population<Gen, Phen>::SelectPhenotypeStrategies(Selector& selector,
                                                 int n) const {
  return gen_eval_->GetPhenotypeStrategies(
      GenotypePopulation<Gen>::SelectGenotypes(selector, n));
}

template <class Gen, class Phen>
std::vector<PhenotypeStrategy<Phen>>
Population<Gen, Phen>::GetAllPhenotypeStrategies() const {
  std::vector<PhenotypeStrategy<Phen>> strats;
  strats.reserve(genes_.size());
  for (const auto& gene : genes_) {
    strats.emplace_back(gen_eval_->GetPhenotypeStrategy(gene));
  }
  return strats;
}

}  // namespace genericga
#endif  // _GENERICGA_POPULATION_H_
