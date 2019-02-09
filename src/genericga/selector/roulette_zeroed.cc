#include "genericga/selector/roulette_zeroed.h"

#include <algorithm>
#include <vector>

namespace genericga {
namespace selector {

std::vector<float> RouletteZeroed::CalculateWeights(
    const std::vector<float>& fitnesses) const {
  auto weights = fitnesses;
  auto min = *std::min_element(std::begin(weights), std::end(weights));
  for (auto& weight : weights) {
    weight -= min;
  }
  return weights;
}

}  // namespace selector
}  // namespace genericga
