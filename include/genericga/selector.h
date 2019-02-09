#ifndef _GENERICGA_SELECTOR_H_
#define _GENERICGA_SELECTOR_H_

#include <vector>

namespace genericga {

class Selector {
 public:
  virtual std::vector<int> SelectIndices(const std::vector<float>& fitnesses,
                                         const std::vector<int>& counts,
                                         int n) = 0;
  virtual ~Selector() {}
};

}  // namespace genericga

#endif  // _GENERICGA_SELECTOR_H_
