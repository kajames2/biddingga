#include "genericga/selector/tournament_mixed.h"

#include <algorithm>
#include <random>
#include <vector>

#include "genericga/vector_ops.h"

namespace genericga {
namespace selector {

TournamentMixed::TournamentMixed(float tourn_size)
    : gen_(std::random_device()()),
      base_tourn_size_(static_cast<int>(tourn_size)),
      frac_extra_(tourn_size - static_cast<int>(tourn_size)) {}

TournamentMixed::TournamentMixed(float tourn_size, int seed)
    : gen_(seed),
      base_tourn_size_(static_cast<int>(tourn_size)),
      frac_extra_(tourn_size - static_cast<int>(tourn_size)) {}

std::vector<int> TournamentMixed::SelectIndices(
    const std::vector<float>& fitnesses, const std::vector<int>& counts,
    int n) {
  ind_dist_ = std::discrete_distribution<>(counts.begin(), counts.end());
  std::vector<int> ind_vec(n);
  std::vector<int> ranks = GetRankingsWithTies(fitnesses, MinRank);
  for (int i = 0; i < n; ++i) {
    ind_vec[i] = TournamentRound(ranks, GenerateTournamentIndices(i, n));
  }
  return ind_vec;
}

int TournamentMixed::TournamentRound(const std::vector<int>& ranks,
                                     std::vector<int> indices) {
  return *std::max_element(
      indices.begin(), indices.end(),
      [&ranks](int i, int j) -> bool { return ranks[i] < ranks[j]; });
}

std::vector<int> TournamentMixed::GenerateTournamentIndices(int cur, int n_draws) {
  int tourn_size = base_tourn_size_;
  if (cur <= frac_extra_ * n_draws) {
    ++tourn_size;
  }
  std::vector<int> tourn_vec(tourn_size);
  for (auto& ind : tourn_vec) {
    ind = ind_dist_(gen_);
  }
  return tourn_vec;
}

}  // namespace selector
}  // namespace genericga
