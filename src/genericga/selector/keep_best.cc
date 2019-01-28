#include "genericga/selector/keep_best.h"

#include <algorithm>
#include <vector>

#include "genericga/fitness_collection.h"

namespace genericga {
namespace selector {

std::vector<int> KeepBest::SelectIndices(const FitnessCollection& col, int n) {
  std::vector<int> selected(n);
  auto orderings = col.GetFitnessOrderings();
  auto it = orderings.rbegin();
  for (auto& val : selected) {
    val = *it;
    ++it;
  }
  return selected;
}

}  // namespace selector
}  // namespace genericga
