#ifndef _GENERICGA_SELECTOR_RANKED_WEIGHTED_H_
#define _GENERICGA_SELECTOR_RANKED_WEIGHTED_H_

#include <cassert>
#include <vector>

#include "genericga/selector/roulette.h"

namespace genericga {
namespace selector {

// Weighted ranking scheme, where a weight of 0 means each index is equally
// likely, while 1 has probability for each proportional to the rank.  0<x<1 is
// a linear combination of the two
class RankedWeighted : public Roulette {
 public:
  RankedWeighted(float weight) : Roulette(), weight_(weight) {
    assert(weight >= 0 && weight <= 1);
  }
  explicit RankedWeighted(float weight, int seed)
      : Roulette(seed), weight_(weight) {
    assert(weight >= 0 && weight <= 1);
  }
  std::vector<float> CalculateWeights(
      const std::vector<float>& fitnesses) const override;

 private:
  float weight_;
};

}  // namespace selector
}  // namespace genericga

#endif  // _GENERICGA_SELECTOR_RANKED_WEIGHTED_H_
