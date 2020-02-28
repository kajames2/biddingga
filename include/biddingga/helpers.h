#ifndef BIDDINGGA_HELPERS_H_
#define BIDDINGGA_HELPERS_H_

#include <random>
#include <vector>

#include "genericga/binary/byte_array_genotype.h"
#include "genericga/multipop/abstract_sub_ga.h"
#include "genericga/multipop/ga.h"
#include "genericga/multipop/sub_ga_adapter.h"

namespace biddingga {

// template <class Environment, class Phen>
// genericga::multipop::GA<Environment> MakeMultipopDriver(
//     const std::vector<std::shared_ptr<
//         genericga::multipop::SubGAAdapter<Environment, Phen>>>& gas,
//     Environment env) {
//   std::vector<std::shared_ptr<genericga::multipop::AbstractSubGA<Environment>>>
//       sub_gas(gas.begin(), gas.end());
//   return genericga::multipop::GA<Environment>(sub_gas, env);
// }

std::vector<unsigned char> RandomByteArray(int n_bytes);

std::vector<genericga::binary::ByteArrayGenotype> RandomGenes(int pop_size,
                                                              int n_bytes);

}  // namespace biddingga

#endif  // BIDDINGGA_HELPERS_H_
