#ifndef GENERICGA_SELECTOR_KEEP_COMMONEST_H_
#define GENERICGA_SELECTOR_KEEP_COMMONEST_H_

#include <vector>

#include "genericga/selector.h"

namespace genericga {
namespace selector {

class KeepCommonest : public Selector {
 public:
  KeepCommonest() {}
  std::vector<int> SelectIndices(const std::vector<float>& fitnesses,
                                 const std::vector<int>& counts,
                                 int n) override;
};

}  // namespace selector
}  // namespace genericga

#endif  // GENERICGA_SELECTOR_KEEP_COMMONEST_H_
