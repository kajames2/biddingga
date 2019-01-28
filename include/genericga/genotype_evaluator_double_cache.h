#ifndef _GENERICGA_GENOTYPE_EVALUATOR_DOUBLE_CACHE_H_
#define _GENERICGA_GENOTYPE_EVALUATOR_DOUBLE_CACHE_H_

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

// GenotypeEvaluatorDoubleCache takes genotypes and stores them as intermediate
// phenotypes. The relationship between a genotype and a phenotype is fixed
// across time, the same genotype will always generate the same phentype for the
// lifetime of the evaluator. The phenotype is converted to a fitness which can
// be updated any time. SetCache sets the pool of genotypes that are stored in
// the evaluator.
template <class Gen, class Phen>
class GenotypeEvaluatorDoubleCache : public AbstractGenotypeEvaluatorCache<Gen, Phen> {
 public:
  explicit GenotypeEvaluatorDoubleCache(GenotypeEvaluator<Gen, Phen> gene_eval);
  GenotypeEvaluatorDoubleCache(std::function<Phen(const Gen&)> phen_conv,
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
  void SetPhenotypeCache(const std::vector<Gen>& genes);
  void SetFitnessCache();
  Phen GetPhenotype(const Gen& gene) const;
  float GetPhenotypeFitness(const Phen& gene) const;

  GenotypeEvaluator<Gen, Phen> gene_eval_;
  mutable std::unordered_map<Gen, Phen, boost::hash<Gen>> phen_cache_;
  mutable std::unordered_map<Phen, float, boost::hash<Phen>> fit_cache_;
};

template <class Gen, class Phen>
GenotypeEvaluatorDoubleCache<Gen, Phen>::GenotypeEvaluatorDoubleCache(
    GenotypeEvaluator<Gen, Phen> gene_eval)
    : gene_eval_(std::move(gene_eval)), phen_cache_(), fit_cache_() {}

template <class Gen, class Phen>
GenotypeEvaluatorDoubleCache<Gen, Phen>::GenotypeEvaluatorDoubleCache(
    std::function<Phen(const Gen&)> phen_conv,
    std::function<float(const Phen&)> fit_calc)
    : gene_eval_(phen_conv, fit_calc) {}

// Create a new cache of strategies, taking already-created strategies from the
// current cache.
template <class Gen, class Phen>
void GenotypeEvaluatorDoubleCache<Gen, Phen>::SetCache(
    const std::vector<Gen>& genes) {
  SetPhenotypeCache(genes);
  SetFitnessCache();
}

template <class Gen, class Phen>
void GenotypeEvaluatorDoubleCache<Gen, Phen>::SetPhenotypeCache(
    const std::vector<Gen>& genes) {
  // Reduce input to unique genes.
  std::unordered_set<Gen, boost::hash<Gen>> temp_genes(genes.begin(),
                                                       genes.end());
  std::vector<Gen> unique_genes(temp_genes.begin(), temp_genes.end());

  // Copy over elements already in the cache.
  std::unordered_map<Gen, Phen, boost::hash<Gen>> new_cache;
  for (const auto& gene : KeyIntersection(unique_genes, phen_cache_)) {
    new_cache.emplace(gene, phen_cache_[gene]);
  }

  // Calculate and insert those not in the cache.
  auto new_genes = KeyDifference(unique_genes, phen_cache_);
  std::vector<Phen> new_phens(new_genes.size());

  #pragma omp parallel for
  for (int i = 0; i < new_genes.size(); ++i) {
    new_phens[i] = gene_eval_.GetPhenotype(new_genes[i]);
  }
  for (int i = 0; i < new_genes.size(); ++i) {
    new_cache.emplace(new_genes[i], new_phens[i]);
  }

  phen_cache_ = new_cache;
}

template <class Gen, class Phen>
void GenotypeEvaluatorDoubleCache<Gen, Phen>::SetFitnessCache() {
  std::unordered_set<Phen, boost::hash<Phen>> temp_phens;
  for (const auto& pair : phen_cache_) {
    temp_phens.emplace(pair.second);
  }

  std::vector<Phen> unique_phens(temp_phens.begin(), temp_phens.end());
  
  // Copy over elements already in the cache.
  std::unordered_map<Phen, float, boost::hash<Phen>> new_cache;
  for (const auto& phen : KeyIntersection(unique_phens, fit_cache_)) {
    new_cache.emplace(phen, fit_cache_[phen]);
  }

  // Calculate and insert those not in the cache.
  auto new_phens = KeyDifference(unique_phens, fit_cache_);
  std::vector<float> new_fits(new_phens.size());
  #pragma omp parallel for
  for (int i = 0; i < new_phens.size(); ++i) {
    new_fits[i] = gene_eval_.GetPhenotypeFitness(new_phens[i]);
  }
  for (int i = 0; i < new_phens.size(); ++i) {
    new_cache.emplace(new_phens[i], new_fits[i]);
  }

  fit_cache_ = new_cache;
}

template <class Gen, class Phen>
void GenotypeEvaluatorDoubleCache<Gen, Phen>::SetFitnessCalculator(
    std::function<float(const Phen&)> fit_calc) {
  gene_eval_.SetFitnessCalculator(fit_calc);
  RecalculateFitnesses();
}

template <class Gen, class Phen>
void GenotypeEvaluatorDoubleCache<Gen, Phen>::RecalculateFitnesses() {
  for (auto& strat : fit_cache_) {
    strat.second = gene_eval_.GetPhenotypeFitness(strat.first);
  }
}

template <class Gen, class Phen>
PhenotypeStrategy<Phen> GenotypeEvaluatorDoubleCache<Gen, Phen>::GetPhenotypeStrategy(
    const Gen& gene) const {
  return PhenotypeStrategy<Phen>{GetPhenotype(gene), GetFitness(gene)};
}

template <class Gen, class Phen>
std::vector<PhenotypeStrategy<Phen>>
GenotypeEvaluatorDoubleCache<Gen, Phen>::GetPhenotypeStrategies(
    const std::vector<Gen>& genes) const {
  std::vector<PhenotypeStrategy<Phen>> strats;
  strats.reserve(genes.size());
  for (const auto& gene : genes) {
    strats.push_back(GetPhenotypeStrategy(gene));
  }
  return strats;
}

template <class Gen, class Phen>
float GenotypeEvaluatorDoubleCache<Gen, Phen>::GetFitness(const Gen& gene) const {
  return GetPhenotypeFitness(GetPhenotype(gene));
}

template <class Gen, class Phen>
Phen GenotypeEvaluatorDoubleCache<Gen, Phen>::GetPhenotype(const Gen& gene) const {
  if (phen_cache_.find(gene) == phen_cache_.end()) {
    phen_cache_.emplace(gene, gene_eval_.GetPhenotype(gene));
  }
  return phen_cache_.at(gene);
}

template <class Gen, class Phen>
float GenotypeEvaluatorDoubleCache<Gen, Phen>::GetPhenotypeFitness(
    const Phen& phen) const {
  if (fit_cache_.find(phen) == fit_cache_.end()) {
    fit_cache_.emplace(phen, gene_eval_.GetPhenotypeFitness(phen));
  }
  return fit_cache_.at(phen);
}

template <class Gen, class Phen>
std::vector<float> GenotypeEvaluatorDoubleCache<Gen, Phen>::GetFitnesses(
    const std::vector<Gen>& genes) const {
  std::vector<float> fits;
  fits.reserve(genes.size());
  for (const auto& gene : genes) {
    fits.push_back(GetFitness(gene));
  }
  return fits;
}

}  // namespace genericga
#endif  // _GENERICGA_GENOTYPE_EVALUATOR_DOUBLE_CACHE_H_
