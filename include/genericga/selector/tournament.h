#ifndef GENERICGA_SELECTOR_TOURNAMENT_H_
#define GENERICGA_SELECTOR_TOURNAMENT_H_

#include <random>
#include <vector>

#include "genericga/selector.h"

namespace genericga {
namespace selector {

class Tournament : public Selector {
 public:
  explicit Tournament(int tourn_size);
  Tournament(int tourn_size, int seed);

  std::vector<int> SelectIndices(const std::vector<float>& fitnesses,
                                 const std::vector<int>& counts,
                                 int n) override;
  int TournamentRound(const std::vector<int>& ranks, std::vector<int> entrants);
  std::vector<int> GenerateTournamentIndices();

 private:
  int tourn_size_;
  std::mt19937 gen_;
  std::discrete_distribution<> dist_;
};

}  // namespace selector
}  // namespace genericga

#endif  // GENERICGA_SELECTOR_TOURNAMENT_H_
