#ifndef _GENERICGA_GENOTYPE_EVALUATOR_SINGLE_CACHE_H_
#define _GENERICGA_GENOTYPE_EVALUATOR_SINGLE_CACHE_H_

#include <algorithm>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "genericga/abstract_genotype_evaluator_cache.h"
#include "genericga/genotype_evaluator.h"
#include "genericga/phenotype_strategy.h"
#include "genericga/vector_ops.h"

namespace genericga {

// GenotypeEvaluatorCache takes genotypes and stores them as intermediate
// phenotypes. The relationship between a genotype and a phenotype is fixed
// across time, the same genotype will always generate the same phentype for the
// lifetime of the evaluator. The phenotype is converted to a fitness which can
// be updated any time. SetCache sets the pool of genotypes that are stored in
// the evaluator.
template <class Gen, class Phen>
class GenotypeEvaluatorSingleCache : public AbstractGenotypeEvaluatorCache<Gen, Phen> {
 public:
  explicit GenotypeEvaluatorSingleCache(GenotypeEvaluator<Gen, Phen> gene_eval);
  GenotypeEvaluatorSingleCache(std::function<Phen(const Gen&)> phen_conv,
                         std::function<float(const Phen&)> fit_calc =
                             [](const Phen&) { return -1.0; });

  void SetCache(const std::vector<Gen>& genotypes) override;
  void SetFitnessCalculator(std::function<float(const Phen&)> fit_calc) override;

  std::vector<float> GetFitnesses(const std::vector<Gen>& genes) const override;
  float GetFitness(const Gen& gene) const override;
  PhenotypeStrategy<Phen> GetPhenotypeStrategy(const Gen& gene) const override;
  std::vector<PhenotypeStrategy<Phen>> GetPhenotypeStrategies(
      const std::vector<Gen>& genes) const override;

 private:
  void RecalculateFitnesses();

  GenotypeEvaluator<Gen, Phen> gene_eval_;
  mutable std::unordered_map<Gen, PhenotypeStrategy<Phen>, boost::hash<Gen>>
      cache_;
};

template <class Gen, class Phen>
GenotypeEvaluatorSingleCache<Gen, Phen>::GenotypeEvaluatorSingleCache(
    GenotypeEvaluator<Gen, Phen> gene_eval)
    : gene_eval_(std::move(gene_eval)), cache_() {}

template <class Gen, class Phen>
GenotypeEvaluatorSingleCache<Gen, Phen>::GenotypeEvaluatorSingleCache(
    std::function<Phen(const Gen&)> phen_conv,
    std::function<float(const Phen&)> fit_calc)
    : gene_eval_(phen_conv, fit_calc) {}

// Create a new cache of strategies, taking already-created strategies from the
// current cache.
template <class Gen, class Phen>
void GenotypeEvaluatorSingleCache<Gen, Phen>::SetCache(
    const std::vector<Gen>& genes) {
  // Reduce input to unique genes.
  std::unordered_set<Gen, boost::hash<Gen>> temp_genes(genes.begin(),
                                                       genes.end());
  std::vector<Gen> unique_genes(temp_genes.begin(), temp_genes.end());

  // Copy over elements already in the cache.
  std::unordered_map<Gen, PhenotypeStrategy<Phen>, boost::hash<Gen>> new_cache;
  for (const auto& gene : KeyIntersection(unique_genes, cache_)) {
    new_cache.emplace(gene, cache_[gene]);
  }

  // Calculate and insert those not in the cache.
  auto new_genes = KeyDifference(unique_genes, cache_);
  std::vector<PhenotypeStrategy<Phen>> new_strats(new_genes.size());
#pragma omp parallel for
  for (int i = 0; i < new_genes.size(); ++i) {
    new_strats[i] = gene_eval_.GetPhenotypeStrategy(new_genes[i]);
  }
  for (int i = 0; i < new_genes.size(); ++i) {
    new_cache.emplace(new_genes[i], new_strats[i]);
  }

  cache_ = new_cache;
}

template <class Gen, class Phen>
void GenotypeEvaluatorSingleCache<Gen, Phen>::SetFitnessCalculator(
    std::function<float(const Phen&)> fit_calc) {
  gene_eval_.SetFitnessCalculator(fit_calc);
  RecalculateFitnesses();
}

template <class Gen, class Phen>
void GenotypeEvaluatorSingleCache<Gen, Phen>::RecalculateFitnesses() {
  for (auto& strat : cache_) {
    strat.second.fitness = gene_eval_.GetPhenotypeFitness(strat.second.phenotype);
  }
}

template <class Gen, class Phen>
PhenotypeStrategy<Phen> GenotypeEvaluatorSingleCache<Gen, Phen>::GetPhenotypeStrategy(
    const Gen& gene) const {
  if (cache_.find(gene) == cache_.end()) {
    cache_.emplace(gene, gene_eval_.GetPhenotypeStrategy(gene));
  }
  return cache_.at(gene);
}

template <class Gen, class Phen>
std::vector<PhenotypeStrategy<Phen>>
GenotypeEvaluatorSingleCache<Gen, Phen>::GetPhenotypeStrategies(
    const std::vector<Gen>& genes) const {
  std::vector<PhenotypeStrategy<Phen>> strats;
  strats.reserve(genes.size());
  for (const auto& gene : genes) {
    strats.push_back(GetPhenotypeStrategy(gene));
  }
  return strats;
}

template <class Gen, class Phen>
float GenotypeEvaluatorSingleCache<Gen, Phen>::GetFitness(const Gen& gene) const {
  return GetPhenotypeStrategy(gene).fitness;
}

template <class Gen, class Phen>
std::vector<float> GenotypeEvaluatorSingleCache<Gen, Phen>::GetFitnesses(
    const std::vector<Gen>& genes) const {
  std::vector<float> fits;
  fits.reserve(genes.size());
  for (const auto& gene : genes) {
    fits.push_back(GetFitness(gene));
  }
  return fits;
}

}  // namespace genericga
#endif  // _GENERICGA_GENOTYPE_EVALUATOR_SINGLE_CACHE_H_
