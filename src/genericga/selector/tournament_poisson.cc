#include "genericga/selector/tournament_poisson.h"

#include <algorithm>
#include <random>
#include <vector>

#include "genericga/vector_ops.h"

namespace genericga {
namespace selector {

TournamentPoisson::TournamentPoisson(float avg_tourn_size)
    : gen_(std::random_device()()),
      extra_size_dist_(avg_tourn_size - 1) {}

TournamentPoisson::TournamentPoisson(float avg_tourn_size, int seed)
    : gen_(seed),
      extra_size_dist_(avg_tourn_size - 1) {}

std::vector<int> TournamentPoisson::SelectIndices(
    const std::vector<float>& fitnesses, const std::vector<int>& counts,
    int n) {
  ind_dist_ = std::discrete_distribution<>(counts.begin(), counts.end());
  std::vector<int> ind_vec(n);
  std::vector<int> ranks = GetRankingsWithTies(fitnesses, MinRank);
  for (int i = 0; i < n; ++i) {
    ind_vec[i] = TournamentRound(ranks, GenerateTournamentIndices());
  }
  return ind_vec;
}

int TournamentPoisson::TournamentRound(const std::vector<int>& ranks,
                                     std::vector<int> indices) {
  return *std::max_element(
      indices.begin(), indices.end(),
      [&ranks](int i, int j) -> bool { return ranks[i] < ranks[j]; });
}

std::vector<int> TournamentPoisson::GenerateTournamentIndices() {
  int n_opponents = extra_size_dist_(gen_);
  std::vector<int> tourn_vec(n_opponents + 1);
  for (auto& ind : tourn_vec) {
    ind = ind_dist_(gen_);
  }
  return tourn_vec;
}

}  // namespace selector
}  // namespace genericga
