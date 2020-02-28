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
#include "genericga/selector/tournament_mixed.h"
#include "numericaldists/interval.h"

namespace biddingga {

struct Configuration1D {
  int id = -1;
  numericaldists::Interval x_range = {0, 1};
  numericaldists::Interval y_range = {0, 1};
  int nx_composites = 10;
  int n_strategies = 100;
  int n_children = 100;
  int nx_segments = 60;
  int bit_precision = 32;
};

struct Configuration2D : public Configuration1D {
  numericaldists::Interval z_range = {0, 1};
  int ny_composites = 10;
  int ny_segments = 60;
};

}  // namespace biddingga

#endif  // BIDDINGGA_CONFIGURATION_H_
