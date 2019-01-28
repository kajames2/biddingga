#include "genericga/binary/bit_mutator.h"

#include <random>

namespace genericga {
namespace binary {

BitMutator::BitMutator(int exp_bit_muts)
    : BitMutator(exp_bit_muts, std::random_device()()) {}

BitMutator::BitMutator(int exp_bit_muts, int seed)
    : gen_(seed), muts_dist_(exp_bit_muts), bit_dist_(0, 0) {}

void BitMutator::operator()(ByteArrayGenotype& gene) {
  if (gene.NBits() != n_bits_) {
    n_bits_ = gene.NBits();
    bit_dist_ = std::uniform_int_distribution<>(0, n_bits_ - 1);
  }
  
  int n_muts = muts_dist_(gen_);
  for (int i = 0; i < n_muts; ++i) {
    gene.Flip(bit_dist_(gen_));
  }
}

}  // namespace binary
}  // namespace genericga
