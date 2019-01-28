#ifndef _GENERICGA_SELECTOR_ROULETTE_SIMPLE_H_
#define _GENERICGA_SELECTOR_ROULETTE_SIMPLE_H_

#include <vector>

#include "genericga/fitness_collection.h"
#include "genericga/selector/roulette.h"

namespace genericga {
namespace selector {

class RouletteSimple : public Roulette {
 public:
  RouletteSimple() : Roulette() {}
  explicit RouletteSimple(int seed) : Roulette(seed) {}

  std::vector<float> CalculateWeights(
      const FitnessCollection& col) const override {
    return col.GetFitnesses();
  }
};

}  // namespace selector
}  // namespace genericga

#endif  // _GENERICGA_SELECTOR_ROULETTE_SIMPLE_H_
