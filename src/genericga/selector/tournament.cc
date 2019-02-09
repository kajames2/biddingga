#include "genericga/selector/tournament.h"

#include <algorithm>
#include <random>
#include <vector>

#include "genericga/vector_ops.h"

namespace genericga {
namespace selector {

Tournament::Tournament(int tourn_size)
    : tourn_size_(tourn_size), gen_(std::random_device()()) {}

Tournament::Tournament(int tourn_size, int seed)
    : tourn_size_(tourn_size), gen_(seed) {}

std::vector<int> Tournament::SelectIndices(const std::vector<float>& fitnesses,
                                           const std::vector<int>& counts,
                                           int n) {
  dist_ = std::discrete_distribution<>(counts.begin(), counts.end());
  std::vector<int> ind_vec(n);
  std::vector<int> ranks = GetRankingsWithTies(fitnesses, MinRank);
  for (int i = 0; i < n; ++i) {
    ind_vec[i] = TournamentRound(ranks, GenerateTournamentIndices());
  }
  return ind_vec;
}

int Tournament::TournamentRound(const std::vector<int>& ranks,
                                std::vector<int> indices) {
  return *std::max_element(
      indices.begin(), indices.end(),
      [&ranks](int i, int j) -> bool { return ranks[i] < ranks[j]; });
}

std::vector<int> Tournament::GenerateTournamentIndices() {
  std::vector<int> tourn_vec(tourn_size_);
  for (auto& ind : tourn_vec) {
    ind = dist_(gen_);
  }
  return tourn_vec;
}

}  // namespace selector
}  // namespace genericga
