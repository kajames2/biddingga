#ifndef GENERICGA_SELECTOR_TOURNAMENT_POISSON_H_
#define GENERICGA_SELECTOR_TOURNAMENT_POISSON_H_

#include <random>
#include <vector>

#include "genericga/selector.h"

namespace genericga {
namespace selector {

class TournamentPoisson : public Selector {
 public:
  explicit TournamentPoisson(float avg_tourn_size);
  TournamentPoisson(float avg_tourn_size, int seed);

  std::vector<int> SelectIndices(const std::vector<float>& fitnesses,
                                 const std::vector<int>& counts,
                                 int n) override;
  int TournamentRound(const std::vector<int>& ranks, std::vector<int> entrants);
  std::vector<int> GenerateTournamentIndices();

 private:
  std::mt19937 gen_;
  std::discrete_distribution<> ind_dist_;
  std::poisson_distribution<> extra_size_dist_;
};

}  // namespace selector
}  // namespace genericga

#endif  // GENERICGA_SELECTOR_TOURNAMENT_POISSON_H_
