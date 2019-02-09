#ifndef _GENERICGA_SELECTOR_ROULETTE_SIMPLE_H_
#define _GENERICGA_SELECTOR_ROULETTE_SIMPLE_H_

#include <vector>

#include "genericga/selector/roulette.h"

namespace genericga {
namespace selector {

class RouletteSimple : public Roulette {
 public:
  RouletteSimple() : Roulette() {}
  explicit RouletteSimple(int seed) : Roulette(seed) {}

  std::vector<float> CalculateWeights(
      const std::vector<float>& fitnesses) const override {
    return fitnesses;
  }
};

}  // namespace selector
}  // namespace genericga

#endif  // _GENERICGA_SELECTOR_ROULETTE_SIMPLE_H_
