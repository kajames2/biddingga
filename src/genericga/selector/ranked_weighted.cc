#include "genericga/selector/ranked_weighted.h"

#include <vector>

#include "genericga/fitness_collection.h"

namespace genericga {
namespace selector {

std::vector<float> RankedWeighted::CalculateWeights(
    const FitnessCollection& col) const {
  auto ranks = col.GetFitnessRankings();
  int size = ranks.size();
  std::vector<float> out_vec(size);
  for (int i = 0; i < size; ++i) {
    out_vec[i] = ((1 - weight_) / size +
                  (2 * ranks[i] * weight_) / (size * (size - 1)));
  }
  return out_vec;
}

}  // namespace selector
}  // namespace genericga
