#ifndef BIDDINGGA_CONFIGURATION_H_
#define BIDDINGGA_CONFIGURATION_H_

#include <memory>

#include "genericga/binary/bit_mutator.h"
#include "genericga/binary/byte_array_genotype.h"
#include "genericga/binary/single_point_crossover.h"
#include "genericga/crossover.h"
#include "genericga/mutator.h"
#include "genericga/selector.h"
#include "genericga/selector/elitism_decorator.h"
#include "genericga/selector/ranked_weighted.h"
#include "numericaldists/interval.h"

namespace biddingga {

struct Configuration1D {
  int id = -1;
  numericaldists::Interval x_range = {0, 1};
  numericaldists::Interval y_range = {0, 1};
  int nx_composites = 5;
  int n_strategies = 100;
  int n_children = 100;
  int nx_segments = 100;
  int bit_precision = 32;
  std::unique_ptr<genericga::Crossover<genericga::binary::ByteArrayGenotype>>
      crossover = std::make_unique<genericga::binary::SinglePointCrossover>();
  std::unique_ptr<genericga::Mutator<genericga::binary::ByteArrayGenotype>>
      mutation = std::make_unique<genericga::binary::BitMutator>(2);
  std::unique_ptr<genericga::Selector> parent_selection =
      std::make_unique<selector::TournamentMixed>(1);
  std::unique_ptr<genericga::Selector> survival =
      std::make_unique<genericga::selector::ElitismDecorator>(
          std::make_unique<genericga::selector::TournamentMixed>(2.2), 2);
};

struct GridConfiguration2D : private Configuration1D {
  Interval z_range = {0, 1};
  int ny_composites = 10;
  int ny_segments = 60;
};

}  // namespace biddingga

#endif  // BIDDINGGA_CONFIGURATION_H_
