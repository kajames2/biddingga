#ifndef _GENERICGA_GENOTYPE_POPULATION_H_
#define _GENERICGA_GENOTYPE_POPULATION_H_

#include "genericga/fitness_collection.h"
#include "genericga/selector.h"

namespace genericga {

template <class Gen>
class GenotypePopulation : public FitnessCollection {
 public:
  virtual void AddGenotypes(std::vector<Gen> genotypes) = 0;
  virtual void SetGenotypes(std::vector<Gen> genotypes) = 0;
  virtual std::vector<Gen> GetAllGenotypes() const = 0;
  Gen SelectGenotype(Selector& selector) const;
  std::vector<Gen> SelectGenotypes(Selector& selector, int n) const;
 protected:
  virtual Gen GetGenotype(int i) const = 0;
  std::vector<Gen> GetGenotypes(std::vector<int> indices) const;
};

template <class Gen>
std::vector<Gen> GenotypePopulation<Gen>::SelectGenotypes(Selector& selector,
                                                          int n) const {
  return GetGenotypes(SelectIndices(selector, n));
}

template <class Gen>
Gen GenotypePopulation<Gen>::SelectGenotype(Selector& selector) const {
  return GetGenotype(SelectIndex(selector));
}

template <class Gen>
std::vector<Gen> GenotypePopulation<Gen>::GetGenotypes(
    std::vector<int> indices) const {
  std::vector<Gen> selected_genes;
  selected_genes.reserve(indices.size());
  for (auto index : indices) {
    selected_genes.emplace_back(GetGenotype(index));
  }
  return selected_genes;
}

}  // namespace genericga

#endif  // _GENERICGA_GENOTYPE_POPULATION_H_
