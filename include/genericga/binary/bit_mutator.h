#ifndef GENERICGA_BINARY_BIT_MUTATOR_H_
#define GENERICGA_BINARY_BIT_MUTATOR_H_

#include <random>

#include "genericga/binary/byte_array_genotype.h"
#include "genericga/mutator.h"

namespace genericga {
namespace binary {

class BitMutator : public Mutator<ByteArrayGenotype> {
 public:
  explicit BitMutator(float exp_bit_muts);
  BitMutator(float exp_bit_muts, int seed);
  void operator()(ByteArrayGenotype& gene) override;

 private:
  std::mt19937 gen_;
  std::poisson_distribution<> muts_dist_;
  std::uniform_int_distribution<> bit_dist_;
  int n_bits_ = 0;
};

}  // namespace binary
}  // namespace genericga

#endif  // GENERICGA_BINARY_BIT_MUTATOR_H_
