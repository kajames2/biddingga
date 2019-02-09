#include "genericga/selector/roulette.h"

#include <random>
#include <vector>

namespace genericga {
namespace selector {

Roulette::Roulette() : gen_(std::random_device()()) {}
Roulette::Roulette(int seed) : gen_(seed) {}

std::vector<int> Roulette::SelectIndices(const std::vector<float>& fitnesses,
                                         const std::vector<int>& counts,
                                         int n) {
  auto weights = CalculateWeights(fitnesses);
  for (int i = 0; i < weights.size(); ++i) {
    weights[i] *= counts[i];
  }
  std::discrete_distribution<> dist(weights.begin(), weights.end());
  std::vector<int> selected;
  for (int i = 0; i < n; ++i) {
    selected.push_back(dist(gen_));
  }
  return selected;
}

}  // namespace selector
}  // namespace genericga
