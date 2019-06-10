#ifndef BIDDINGGA_FUNCTION_GA_CONFIGURATION_H_
#define BIDDINGGA_FUNCTION_GA_CONFIGURATION_H_

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

struct FunctionGAConfiguration {
  numericaldists::Interval value_range = {0, 1};
  numericaldists::Interval bid_range = {0, 1};
  int n_strategies = 1000;
  int n_children = 1000;
  int n_segments = 30;
  int bit_precision = 32;
  std::unique_ptr<genericga::Crossover<genericga::binary::ByteArrayGenotype>>
      crossover = std::make_unique<genericga::binary::SinglePointCrossover>();
  std::unique_ptr<genericga::Mutator<genericga::binary::ByteArrayGenotype>>
      mutation = std::make_unique<genericga::binary::BitMutator>(3);
  std::unique_ptr<genericga::Selector> survival =
      std::make_unique<genericga::selector::ElitismDecorator>(
          std::make_unique<genericga::selector::RankedWeighted>(0.6), 5);
};

}  // namespace biddingga

#endif  // BIDDINGGA_FUNCTION_GA_CONFIGURATION_H_
