#ifndef GENERICGA_POPULATION_H_
#define GENERICGA_POPULATION_H_

#include <algorithm>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>
#include <utility>
#include <vector>

#include "genericga/genotype_population.h"
#include "genericga/phenotype_strategy.h"

namespace genericga {

template <class Gen, class Phen>
class Population : public GenotypePopulation<Gen> {
 public:
  Population(
      std::function<Phen(const Gen&)> phen_conv,
      std::function<std::vector<float>(const std::vector<Phen>&)> fit_calc,
      std::vector<Gen> genes);

  void SetFitnessCalculator(
      std::function<std::vector<float>(const std::vector<Phen>&)> fit_calc);

  void AddGenotypes(std::vector<Gen> genes) override;
  void Survival(Selector& selector, int n) override;

  Gen SelectGenotype(Selector& selector) const override;
  std::vector<Gen> SelectGenotypes(Selector& selector, int n) const override;
  PhenotypeStrategy<Phen> SelectPhenotypeStrategy(Selector& selector) const;
  std::vector<PhenotypeStrategy<Phen>> SelectPhenotypeStrategies(
      Selector& selector, int n) const;
  int GetTotalCount() const {
    return std::accumulate(counts_.begin(), counts_.end(), 0);
  }

  std::vector<Gen> GetGenotypes() const override { return genes_; }
  std::vector<Phen> GetPhenotypes() const { return phens_; }
  std::vector<float> GetFitnesses() const override { return fits_; }
  std::vector<int> GetCounts() const { return counts_; }
  std::vector<PhenotypeStrategy<Phen>> GetPhenotypeStrategies() const;

 private:
  void RemoveDead();

  std::vector<Gen> genes_;
  std::vector<Phen> phens_;
  std::vector<float> fits_;
  std::vector<int> counts_;
  std::function<Phen(const Gen&)> phen_conv_;
  std::function<std::vector<float>(const std::vector<Phen>&)> fit_calc_;
};

template <class Gen, class Phen>
Population<Gen, Phen>::Population(
    std::function<Phen(const Gen&)> phen_conv,
    std::function<std::vector<float>(const std::vector<Phen>&)> fit_calc,
    std::vector<Gen> genes)
    : genes_(std::move(genes)), phen_conv_(std::move(phen_conv)) {
  phens_.reserve(genes_.size());
  counts_ = std::vector<int>(genes_.size(), 1);
  for (const auto& gene : genes_) {
    phens_.push_back(phen_conv_(gene));
  }
  SetFitnessCalculator(fit_calc);
}

template <class Gen, class Phen>
void Population<Gen, Phen>::SetFitnessCalculator(
    std::function<std::vector<float>(const std::vector<Phen>&)> fit_calc) {
  fit_calc_ = std::move(fit_calc);
  fits_ = fit_calc_(phens_);
}

template <class Gen, class Phen>
void Population<Gen, Phen>::AddGenotypes(std::vector<Gen> new_genes) {
  std::vector<Phen> new_phens;
  new_phens.reserve(new_genes.size());
  std::transform(new_genes.begin(), new_genes.end(),
                 std::back_inserter(new_phens), phen_conv_);
  auto new_fits = fit_calc_(new_phens);

  std::copy(std::make_move_iterator(new_genes.begin()),
            std::make_move_iterator(new_genes.end()),
            std::back_inserter(genes_));
  std::copy(std::make_move_iterator(new_phens.begin()),
            std::make_move_iterator(new_phens.end()),
            std::back_inserter(phens_));
  std::copy(new_fits.begin(), new_fits.end(), std::back_inserter(fits_));

  std::fill_n(std::back_inserter(counts_), new_genes.size(), 1);
}

template <class Gen, class Phen>
void Population<Gen, Phen>::Survival(Selector& selector, int n) {
  auto inds = selector.SelectIndices(fits_, counts_, n);
  std::fill(counts_.begin(), counts_.end(), 0);
  for (int ind : inds) {
    ++counts_[ind];
  }
  RemoveDead();
}

// 0-count strategies are "Dead".  Move non-dead from the end to take their
// place, and then truncate the vectors
template <class Gen, class Phen>
void Population<Gen, Phen>::RemoveDead() {
  int size = counts_.size();
  for (int i = 0; i < size; ++i) {
    // Find an open spot due to a dead strategy
    if (counts_[i] == 0) {
      // Find the first non-dead strategy from the end
      while (i != size - 1 && counts_[size - 1] == 0) {
        --size;
      }
      // Fill in the hole with the non-dead strategy
      if (i != size - 1) {
        genes_[i] = std::move(genes_[size - 1]);
        phens_[i] = std::move(phens_[size - 1]);
        fits_[i] = fits_[size - 1];
        counts_[i] = counts_[size - 1];
        --size;
      }
    }
  }
  // Truncate the lists
  genes_.erase(genes_.begin() + size, genes_.end());
  phens_.erase(phens_.begin() + size, phens_.end());
  fits_.resize(size);
  counts_.resize(size);
}

template <class Gen, class Phen>
Gen Population<Gen, Phen>::SelectGenotype(Selector& selector) const {
  return SelectGenotypes(selector, 1)[0];
}

template <class Gen, class Phen>
std::vector<Gen> Population<Gen, Phen>::SelectGenotypes(Selector& selector,
                                                        int n) const {
  auto inds = selector.SelectIndices(fits_, counts_, n);
  std::vector<Gen> out_genes;
  out_genes.reserve(n);
  for (int i : inds) {
    out_genes.push_back(genes_[i]);
  }
  return genes_;
}

template <class Gen, class Phen>
PhenotypeStrategy<Phen> Population<Gen, Phen>::SelectPhenotypeStrategy(
    Selector& selector) const {
  return SelectPhenotypeStrategies(selector, 1)[0];
}

template <class Gen, class Phen>
std::vector<PhenotypeStrategy<Phen>>
Population<Gen, Phen>::SelectPhenotypeStrategies(Selector& selector,
                                                 int n) const {
  auto inds = selector.SelectIndices(fits_, counts_, n);
  std::vector<PhenotypeStrategy<Phen>> strats;
  strats.reserve(n);
  for (int ind : inds) {
    strats.push_back({phens_[ind], fits_[ind]});
  }
  return strats;
}

template <class Gen, class Phen>
std::vector<PhenotypeStrategy<Phen>>
Population<Gen, Phen>::GetPhenotypeStrategies() const {
  std::vector<PhenotypeStrategy<Phen>> strats;
  strats.reserve(phens_.size());
  for (int i = 0; i < phens_.size(); ++i) {
    strats.push_back({phens_[i], fits_[i]});
  }
  return strats;
}

}  // namespace genericga
#endif  // GENERICGA_POPULATION_H_
