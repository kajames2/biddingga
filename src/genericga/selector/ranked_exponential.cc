#include "genericga/selector/ranked_exponential.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "genericga/fitness_collection.h"

namespace genericga{
namespace selector{

std::vector<float> RankedExponential::CalculateWeights(
    const FitnessCollection& col) const {
  auto ranks = col.GetFitnessRankings();
  std::transform(ranks.begin(), ranks.end(), ranks.begin(),
                 [](float rank) -> float { return 1 - std::exp(-rank); });
  return ranks;
}

}  // namespace genericga
}  // namespace selector
