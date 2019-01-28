#ifndef _GENERICGA_SELECTOR_TOURNAMENT_H_
#define _GENERICGA_SELECTOR_TOURNAMENT_H_

#include <random>
#include <vector>

#include "genericga/fitness_collection.h"
#include "genericga/selector.h"

namespace genericga {
namespace selector {

class Tournament : public Selector {
 public:
  Tournament(int tourn_size);
  Tournament(int tourn_size, int seed);

  std::vector<int> SelectIndices(const FitnessCollection& col, int n) override;
  int TournamentRound(const std::vector<float>& col, std::vector<int> entrants);
  std::vector<int> GenerateTournamentIndices();

 private:
  int tourn_size_;
  std::mt19937 gen_;
  std::uniform_int_distribution<> dist;
};

}  // namespace selector
}  // namespace genericga

#endif  // _GENERICGA_SELECTOR_TOURNAMENT_H_
