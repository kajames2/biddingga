#include "genericga/selector/keep_commonest.h"

#include <algorithm>
#include <vector>

#include "genericga/vector_ops.h"

namespace genericga {
namespace selector {

std::vector<int> KeepCommonest::SelectIndices(const std::vector<float>& fitnesses,
                                         const std::vector<int>& counts,
                                         int n) {
  std::vector<int> selected(n);
  auto orderings = GetOrderings(counts);
  auto it = orderings.rbegin();
  int count = counts[*it];
  for (auto& val : selected) {
    val = *it;
    --count;
    if (count == 0) {
      ++it;
      count = counts[*it];
    }
  }
  return selected;
}

}  // namespace selector
}  // namespace genericga
