#include "genericga/selector/ranked_weighted.h"

#include <vector>

#include "genericga/vector_ops.h"

namespace genericga {
namespace selector {

std::vector<float> RankedWeighted::CalculateWeights(
    const std::vector<float>& fitnesses) const {
  auto ranks = GetRankingsWithTies(fitnesses, AverageRank);
  int size = ranks.size();
  std::vector<float> out_vec(size);
  for (int i = 0; i < size; ++i) {
    out_vec[i] = ((1 - weight_) * 1.0 / size) +
                 (weight_ * (2 * ranks[i]) / (size * (size - 1)));
  }
  return out_vec;
}

}  // namespace selector
}  // namespace genericga
