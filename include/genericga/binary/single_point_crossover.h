#ifndef _GENERICGA_BINARY_SINGLE_POINT_CROSSOVER_H_
#define _GENERICGA_BINARY_SINGLE_POINT_CROSSOVER_H_

#include <climits>
#include <iostream>
#include <random>

#include "genericga/binary/byte_array_genotype.h"
#include "genericga/crossover.h"

namespace genericga {
namespace binary {

class SinglePointCrossover : public Crossover<ByteArrayGenotype> {
 public:
  SinglePointCrossover();
  explicit SinglePointCrossover(int seed);
  void operator()(ByteArrayGenotype& gene1, ByteArrayGenotype& gene2) override;

 private:
  std::mt19937 gen_;
  int n_bits_ = 0;
  std::uniform_int_distribution<> dist;
};

void CrossoverGenes(ByteArrayGenotype& gene1, ByteArrayGenotype& gene2, int bit1,
               int bit2);

}  // namespace binary
}  // namespace genericga

#endif  // _GENERICGA_BINARY_SINGLE_POINT_CROSSOVER_H_
