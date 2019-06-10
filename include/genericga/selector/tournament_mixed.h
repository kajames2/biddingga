#ifndef GENERICGA_SELECTOR_TOURNAMENT_MIXED_H_
#define GENERICGA_SELECTOR_TOURNAMENT_MIXED_H_

#include <random>
#include <vector>

#include "genericga/selector.h"

namespace genericga {
namespace selector {

class TournamentMixed : public Selector {
 public:
  explicit TournamentMixed(float tourn_size);
  TournamentMixed(float tourn_size, int seed);

  std::vector<int> SelectIndices(const std::vector<float>& fitnesses,
                                 const std::vector<int>& counts,
                                 int n) override;
  int TournamentRound(const std::vector<int>& ranks, std::vector<int> entrants);
  std::vector<int> GenerateTournamentIndices(int cur, int n_draws);

 private:
  int base_tourn_size_;
  float frac_extra_;
  std::mt19937 gen_;
  std::discrete_distribution<> ind_dist_;
};

}  // namespace selector
}  // namespace genericga

#endif  // GENERICGA_SELECTOR_TOURNAMENT_MIXED_H_
