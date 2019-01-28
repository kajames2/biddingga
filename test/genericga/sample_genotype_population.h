#ifndef _SAMPLE_GENOTYPE_POPULATION_H_
#define _SAMPLE_GENOTYPE_POPULATION_H_

#include <vector>

#include "genericga/genotype_population.h"

namespace gatests {
class SampleGenotypePopulation : public genericga::GenotypePopulation<int> {
 public:
  SampleGenotypePopulation()
      : fits(std::vector<float>{1, 3, 5, -2}), genes(std::vector<int>{6, 7, 8, 2}) {}
  std::vector<float> GetFitnesses() const override { return fits; }
  float GetFitness(int i) const override { return fits[i]; }

  void AddGenotypes(std::vector<int> genotypes) override {}
  void SetGenotypes(std::vector<int> genotypes) override {}
  std::vector<int> GetAllGenotypes() const override {return genes;}

 protected:
  int GetGenotype(int i) const override {return genes[i];}
  
 private:
  std::vector<float> fits;
  std::vector<int> genes;
};
}  // namespace gatests

#endif  // _SAMPLE_GENOTYPE_POPULATION_H_
