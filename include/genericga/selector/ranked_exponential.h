#ifndef _GENERICGA_SELECTOR_RANKED_EXPONENTIAL_H_
#define _GENERICGA_SELECTOR_RANKED_EXPONENTIAL_H_

#include <vector>

#include "genericga/fitness_collection.h"
#include "genericga/selector/roulette.h"

namespace genericga {
namespace selector {

class RankedExponential : public Roulette {
 public:
  RankedExponential() : Roulette() {}
  explicit RankedExponential(int seed) : Roulette(seed) {}
  std::vector<float> CalculateWeights(
      const FitnessCollection& col) const override;
};

}  // namespace selector
}  // namespace genericga

#endif  // _GENERICGA_SELECTOR_RANKED_EXPONENTIAL_H_
