#ifndef GENERICGA_SELECTOR_RANKED_EXPONENTIAL_H_
#define GENERICGA_SELECTOR_RANKED_EXPONENTIAL_H_

#include <vector>

#include "genericga/selector/roulette.h"

namespace genericga {
namespace selector {

class RankedExponential : public Roulette {
 public:
  RankedExponential() : Roulette() {}
  explicit RankedExponential(int seed) : Roulette(seed) {}
  std::vector<float> CalculateWeights(
      const std::vector<float>& fitnesses) const override;
};

}  // namespace selector
}  // namespace genericga

#endif  // GENERICGA_SELECTOR_RANKED_EXPONENTIAL_H_
