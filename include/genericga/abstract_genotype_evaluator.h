#ifndef _GENERICGA_ABSTRACT_GENOTYPE_EVALUATOR_H_
#define _GENERICGA_ABSTRACT_GENOTYPE_EVALUATOR_H_

#include <vector>

namespace genericga {

template <class Gen>
class AbstractGenotypeEvaluator {
 public:
  virtual std::vector<float> GetFitnesses(const std::vector<Gen>& genes) const = 0;
  virtual float GetFitness(const Gen& gene) const = 0;
  virtual ~AbstractGenotypeEvaluator() {}
};

}  // namespace genericga
#endif  // _GENERICGA_ABSTRACT_GENOTYPE_EVALUATOR_H_
