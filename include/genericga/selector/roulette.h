#ifndef _GENERICGA_SELECTOR_ROULETTE_H_
#define _GENERICGA_SELECTOR_ROULETTE_H_

#include <random>
#include <vector>

#include "genericga/fitness_collection.h"
#include "genericga/selector.h"

namespace genericga {
namespace selector {

class Roulette : public Selector {
 public:
  Roulette();
  explicit Roulette(int seed);
  std::vector<int> SelectIndices(const FitnessCollection& col, int n) override;
  virtual std::vector<float> CalculateWeights(
      const FitnessCollection& col) const = 0;

 private:
  std::mt19937 gen_;
};

}  // namespace selector
}  // namespace genericga

#endif  // _GENERICGA_SELECTOR_ROULETTE_H_
