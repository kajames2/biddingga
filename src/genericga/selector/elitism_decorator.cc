#include "genericga/selector/elitism_decorator.h"

#include <memory>
#include <vector>

#include "genericga/fitness_collection.h"
#include "genericga/selector.h"
#include "genericga/selector/keep_best.h"

namespace genericga {
namespace selector {

ElitismDecorator::ElitismDecorator(std::unique_ptr<Selector> sel, int n_elites)
    : sel_(std::move(sel)),
      elite_sel_(std::make_unique<KeepBest>()),
      n_elites_(n_elites) {}

std::vector<int> ElitismDecorator::SelectIndices(const FitnessCollection& col,
                                                 int n) {
  int actual_n_elites = std::min(n, n_elites_);
  auto elites = elite_sel_->SelectIndices(col, actual_n_elites);
  auto rest = sel_->SelectIndices(col, n - actual_n_elites);
  rest.insert(rest.end(), elites.begin(), elites.end());
  return rest;
}

}  // namespace selector
}  // namespace genericga
