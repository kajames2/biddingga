#ifndef _GENERICGA_SELECTOR_ROULETTE_ZEROED_H_
#define _GENERICGA_SELECTOR_ROULETTE_ZEROED_H_

#include <vector>

#include "genericga/selector/roulette.h"

namespace genericga {
namespace selector {

class RouletteZeroed : public Roulette {
 public:
  RouletteZeroed() : Roulette() {}
  explicit RouletteZeroed(int seed) : Roulette(seed) {}

  std::vector<float> CalculateWeights(
      const std::vector<float>& fitnesses) const override;
};

}  // namespace selector
}  // namespace genericga

#endif  // _GENERICGA_ROULETTE_ZEROED_SELECTOR_H_
