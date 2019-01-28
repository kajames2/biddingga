#ifndef _GENERICGA_SELECTOR_H_
#define _GENERICGA_SELECTOR_H_

#include <vector>

namespace genericga {

class FitnessCollection;

class Selector {
 public:
  virtual std::vector<int> SelectIndices(const FitnessCollection& col,
                                         int n) = 0;
  virtual ~Selector() {}
};

}  // namespace genericga

#endif  // _GENERICGA_SELECTOR_H_
