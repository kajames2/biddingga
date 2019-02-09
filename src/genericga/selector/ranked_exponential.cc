#include "genericga/selector/ranked_exponential.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "genericga/vector_ops.h"

namespace genericga{
namespace selector{

std::vector<float> RankedExponential::CalculateWeights(
    const std::vector<float>& fitnesses) const {
  auto ranks = GetRankingsWithTies(fitnesses, AverageRank);
  std::transform(ranks.begin(), ranks.end(), ranks.begin(),
                 [](float rank) -> float { return 1 - std::exp(-rank); });
  return ranks;
}

}  // namespace genericga
}  // namespace selector
