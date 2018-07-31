#include <array>
#include <map>

#include <boost/math/distributions/uniform.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include "genericga/binary/bit_mutator.h"
#include "genericga/binary/byte_array_genotype.h"
#include "genericga/binary/encoding.h"
#include "genericga/binary/single_point_crossover.h"
#include "genericga/multipop/ga.h"
#include "genericga/multipop/single_population_ga_to_sub_ga_adapter.h"
#include "genericga/selector/elitism_decorator.h"
#include "genericga/selector/keep_best.h"
#include "genericga/selector/ranked_weighted.h"
#include "genericga/selector/tournament.h"
#include "genericga/single_population_ga.h"


namespace biddingga {

using namespace genericga;
using namespace bidding;

using Gen = binary::ByteArrayGenotype;
using Phen = PiecewiseLinearFunction;

template <class Phen>
SinglePopulationGA<binary::ByteArrayGenotype, Phen> BinaryGA(
    std::function<Phen(const binary::ByteArrayGenotype&)> phen_conv,
    int pop_size, int n_bits,
    std::function<float(const Phen&)> fit = [](const Phen&) { return -1.0; }) {
  std::vector<Gen> genes;
  genes.reserve(pop_size);
  for (int i = 0; i < pop_size; ++i) {
    genes.emplace_back(
        std::vector<unsigned char>((n_bits + CHAR_BIT - 1) / CHAR_BIT));
  }
  Population<Gen, Phen> init_pop(phen_conv, fit, genes);

  auto children_fact = std::make_unique<ChildrenFactory<Gen>>(
      std::make_unique<binary::SinglePointCrossover>(),
      std::make_unique<binary::BitMutator>(3));

  return SinglePopulationGA<binary::ByteArrayGenotype, Phen>(
      std::move(init_pop), std::move(children_fact));
}

struct PhenConvert {
  PiecewiseLinearFunction operator()(const Gen& gene) const {
    std::vector<float> ys = gene.ToFloatArray(nums_);
    std::sort(ys.begin(), ys.end());
    return PiecewiseLinearFunction(ys, range_);
  }
  std::vector<binary::Encoding> nums_;
  Range range_;
};

int main(int argc, char** argv) {
  std::vector<binary::Encoding> encodings(40, binary::Encoding{32, 0, 150});
  int n_bits = 0;
  for (const auto& en : encodings) {
    n_bits += en.bit_precision;
  }
  PhenConvert conv{encodings, Range{0, 100}};
  std::vector<uniform_distribution<>> dists{uniform_distribution<>(0, 100),
                                            uniform_distribution<>(0, 100),
                                            uniform_distribution<>(50, 150)};
  FirstPriceAuctionEnvironment auction(dists);
  std::vector<
      std::shared_ptr<multipop::AbstractSubGA<FirstPriceAuctionEnvironment>>>
      sub_gas;
  std::vector<std::shared_ptr<multipop::SinglePopulationGAToSubGAAdapter<
      FirstPriceAuctionEnvironment, Phen>>>
      ind_gas;

  int pop_size = 200;
  for (int i = 0; i < dists.size(); ++i) {
    auto ga = std::make_unique<SinglePopulationGA<Gen, Phen>>(
        BinaryGA<Phen>(conv, pop_size, n_bits));

    auto subga = std::make_shared<multipop::SinglePopulationGAToSubGAAdapter<
        FirstPriceAuctionEnvironment, Phen>>(std::move(ga), i);
    sub_gas.push_back(subga);
    ind_gas.push_back(subga);
  }
  multipop::GA<FirstPriceAuctionEnvironment> ga(sub_gas, auction);

  ga.RunRound(200);
  std::vector<Phen> phens;
  for (int i = 0; i < dists.size(); ++i) {
    phens.push_back(ind_gas[i]->GetBestStrategy().phenotype);
  }
  for (int i = 0; i <= 100; ++i) {
    std::cout << i << ": ";
    for (int player = 0; player < phens.size(); ++player) {
      if (pdf(dists[player], i) == 0) {
        std::cout << "---" << '\t';
      } else {
        std::cout << phens[player].GetBid(i) << '\t';
      }
    }
    std::cout << '\n';
  }
  std::cout << std::endl;

  return 0;
}

}  // namespace biddingga
