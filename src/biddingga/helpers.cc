#include "biddingga/helpers.h"
#include <vector>

namespace biddingga {

std::vector<unsigned char> RandomByteArray(int n_bytes) {
  static auto generator = std::mt19937(std::random_device()());
  std::uniform_int_distribution<int> dist(0, UCHAR_MAX);
  std::vector<unsigned char> rand_gene(n_bytes);
  for (int j = 0; j < n_bytes; ++j) {
    rand_gene[j] = static_cast<unsigned char>(dist(generator));
  }
  return rand_gene;
}

std::vector<genericga::binary::ByteArrayGenotype> RandomGenes(int pop_size,
                                                              int n_bytes) {
  std::vector<genericga::binary::ByteArrayGenotype> genes;
  for (int i = 0; i < pop_size; ++i) {
    genes.emplace_back(RandomByteArray(n_bytes));
  }
  return genes;
}

}  // namespace biddingga
