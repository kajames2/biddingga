#ifndef _GENERICGA_GENOTYPE_POPULATION_H_
#define _GENERICGA_GENOTYPE_POPULATION_H_

#include "genericga/selector.h"

namespace genericga {

template <class Gen>
class GenotypePopulation {
 public:
  virtual void AddGenotypes(std::vector<Gen> genotypes) = 0;
  virtual void Survival(Selector& selector, int n) = 0;
  virtual std::vector<Gen> GetGenotypes() const = 0;
  virtual Gen SelectGenotype(Selector& selector) const = 0;
  virtual std::vector<Gen> SelectGenotypes(Selector& selector, int n) const = 0;
  virtual std::vector<float> GetFitnesses() const = 0;
};

}  // namespace genericga

#endif  // _GENERICGA_GENOTYPE_POPULATION_H_
