#ifndef _GENERICGA_SELECTOR_KEEP_BEST_H_
#define _GENERICGA_SELECTOR_KEEP_BEST_H_

#include <vector>

#include "genericga/fitness_collection.h"
#include "genericga/selector.h"

namespace genericga {
namespace selector {

class KeepBest : public Selector {
 public:
  KeepBest() {}
  std::vector<int> SelectIndices(const FitnessCollection& col, int n) override;
};

}  // namespace selector
}  // namespace genericga

#endif  // _GENERICGA_SELECTOR_KEEP_BEST_H_
