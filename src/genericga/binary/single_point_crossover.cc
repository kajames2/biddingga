#include "genericga/binary/single_point_crossover.h"

#include <random>

namespace genericga {
namespace binary {

SinglePointCrossover::SinglePointCrossover()
    : SinglePointCrossover(std::random_device()()) {}

SinglePointCrossover::SinglePointCrossover(int seed) : gen_(seed), dist(0, 0) {}

void SinglePointCrossover::operator()(ByteArrayGenotype& gene1,
                                      ByteArrayGenotype& gene2) {
  if (gene1.NBits() != n_bits_ || gene2.NBits() != n_bits_) {
    n_bits_ = std::min(gene1.NBits(), gene2.NBits());
    dist = std::uniform_int_distribution<>(0, n_bits_ - 1);
  }
  int bit1 = dist(gen_);
  int bit2 = dist(gen_);
  CrossoverGenes(gene1, gene2, bit1, bit2);
}

void CrossoverGenes(ByteArrayGenotype& gene1, ByteArrayGenotype& gene2,
                    int bit1, int bit2) {
  bit1 < bit2 ? SwapBits(gene1, gene2, bit1, bit2 - bit1 + 1)
              : SwapBits(gene1, gene2, bit2, bit1 - bit2 + 1);
}

}  // namespace binary
}  // namespace genericga
