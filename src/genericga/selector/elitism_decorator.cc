#include "genericga/selector/elitism_decorator.h"

#include <memory>
#include <vector>

#include "genericga/selector.h"
#include "genericga/selector/keep_best.h"

namespace genericga {
namespace selector {

ElitismDecorator::ElitismDecorator(std::unique_ptr<Selector> sel, int n_elites)
    : sel_(std::move(sel)),
      elite_sel_(std::make_unique<KeepBest>()),
      n_elites_(n_elites) {}

std::vector<int> ElitismDecorator::SelectIndices(
    const std::vector<float>& fitnesses, const std::vector<int>& counts,
    int n) {
  int actual_n_elites = std::min(n, n_elites_);
  auto elites = elite_sel_->SelectIndices(fitnesses, counts, actual_n_elites);
  auto rest = sel_->SelectIndices(fitnesses, counts, n - actual_n_elites);
  rest.insert(rest.end(), elites.begin(), elites.end());
  return rest;
}

}  // namespace selector
}  // namespace genericga
