#include "genericga/selector/roulette_zeroed.h"

#include <algorithm>
#include <vector>

#include "genericga/fitness_collection.h"

namespace genericga {
namespace selector {

std::vector<float> RouletteZeroed::CalculateWeights(
    const FitnessCollection& col) const {
  auto weights = col.GetFitnesses();
  auto min = *std::min_element(std::begin(weights), std::end(weights));
  std::transform(weights.begin(), weights.end(), weights.begin(),
                 [min](float weight) -> float { return weight - min; });
  return weights;
}

}  // namespace selector
}  // namespace genericga
