#ifndef _GENERICGA_SELECTOR_RANKED_WEIGHTED_H_
#define _GENERICGA_SELECTOR_RANKED_WEIGHTED_H_

#include <vector>

#include "genericga/fitness_collection.h"
#include "genericga/selector/roulette.h"

namespace genericga {
namespace selector {

// Weighted ranking scheme, where a weight of 0 means each index is equally
// likely, while 1 has probability for each proportional to the rank.  0<x<1 is
// a linear combination of the two
class RankedWeighted : public Roulette {
 public:
  RankedWeighted(float weight) : Roulette(), weight_(weight) {}
  explicit RankedWeighted(float weight, int seed)
      : Roulette(seed), weight_(weight) {}
  std::vector<float> CalculateWeights(
      const FitnessCollection& col) const override;

 private:
  float weight_;
};

}  // namespace selector
}  // namespace genericga

#endif  // _GENERICGA_SELECTOR_RANKED_WEIGHTED_H_
