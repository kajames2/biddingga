#ifndef _GENERICGA_SELECTOR_KEEP_BEST_H_
#define _GENERICGA_SELECTOR_KEEP_BEST_H_

#include <vector>

#include "genericga/selector.h"

namespace genericga {
namespace selector {

class KeepBest : public Selector {
 public:
  KeepBest() {}
  std::vector<int> SelectIndices(const std::vector<float>& fitnesses,
                                 const std::vector<int>& counts,
                                 int n) override;
};

}  // namespace selector
}  // namespace genericga

#endif  // _GENERICGA_SELECTOR_KEEP_BEST_H_
