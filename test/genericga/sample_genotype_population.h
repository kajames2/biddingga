#ifndef _SAMPLE_GENOTYPE_POPULATION_H_
#define _SAMPLE_GENOTYPE_POPULATION_H_

#include <vector>

#include "genericga/genotype_population.h"
#include "genericga/selector.h"

namespace gatests {
class SampleGenotypePopulation : public genericga::GenotypePopulation<int> {
 public:
  SampleGenotypePopulation()
      : fits(std::vector<float>{1, 3, 5, -2}),
        genes(std::vector<int>{6, 7, 8, 2}) {}
  void AddGenotypes(std::vector<int> genotypes) override {}
  void Survival(genericga::Selector& selector, int n) override {}
  int SelectGenotype(genericga::Selector& selector) const override {
    return SelectGenotypes(selector, 1)[0];
  }
  std::vector<int> SelectGenotypes(genericga::Selector& selector,
                                   int n) const override {
    auto inds =
        selector.SelectIndices(fits, std::vector<int>(fits.size(), 1), n);
    std::vector<int> out_genes;
    for (int i : inds) {
      out_genes.push_back(genes[i]);
    }
    return out_genes;
  }

  std::vector<float> GetFitnesses() const override { return fits; }
  std::vector<int> GetGenotypes() const override { return genes; }

 protected:
 private:
  std::vector<float> fits;
  std::vector<int> genes;
};
}  // namespace gatests

#endif  // _SAMPLE_GENOTYPE_POPULATION_H_
