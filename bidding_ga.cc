#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/uniform.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include "bidding/all_pay_auction.h"
#include "bidding/distribution.h"
#include "bidding/equal_length_piecewise_function.h"
#include "bidding/first_price_auction.h"
#include "bidding/piecewise_linear_function.h"
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

using namespace genericga;
using namespace bidding;
using namespace boost::math;

using Gen = binary::ByteArrayGenotype;
using Phen = EqualLengthPiecewiseFunction;

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

struct PhenConvertRaw {
  EqualLengthPiecewiseFunction operator()(const Gen& gene) const {
    std::vector<float> ys = gene.ToFloatArray(nums_);
    return EqualLengthPiecewiseFunction(ys, interval_);
  }
  std::vector<binary::Encoding> nums_;
  Interval interval_;
};

struct PhenConvertSort {
  EqualLengthPiecewiseFunction operator()(const Gen& gene) const {
    std::vector<float> ys = gene.ToFloatArray(nums_);
    std::sort(ys.begin(), ys.end());
    return EqualLengthPiecewiseFunction(ys, interval_);
  }
  std::vector<binary::Encoding> nums_;
  Interval interval_;
};

struct PhenConvertFlooring {
  EqualLengthPiecewiseFunction operator()(const Gen& gene) const {
    std::vector<float> ys = gene.ToFloatArray(nums_);
    float floor = ys[0];
    for (auto& val : ys) {
      if (val < floor) {
        val = floor;
      }
      floor = val;
    }
    return EqualLengthPiecewiseFunction(ys, interval_);
  }
  std::vector<binary::Encoding> nums_;
  Interval interval_;
};

struct PhenConvertCumulative {
  EqualLengthPiecewiseFunction operator()(const Gen& gene) const {
    std::vector<float> ys = gene.ToFloatArray(nums_);
    float tot = 0;
    for (auto& val : ys) {
      tot += val;
      if (tot > 1) {
        tot = 1;
       }
      val = tot;
    }
    return EqualLengthPiecewiseFunction(ys, interval_);
  }
  std::vector<binary::Encoding> nums_;
  Interval interval_;
};

void RunFirstPrice() {
  using PhenConvert = PhenConvertSort;
  std::vector<Distribution> dists{
      uniform_distribution<>(0, 100),
      uniform_distribution<>(0, 100),
      uniform_distribution<>(0, 100)};
      // exponential_distribution<>(0.04),
      // normal_distribution<>(30, 4), normal_distribution<>(70, 4),
      // normal_distribution<>(110, 4)};
  float minValue = lower(dists);
  float maxValue = upper(dists);
  std::vector<binary::Encoding> encodings(30,
                                          binary::Encoding{32, 0, maxValue});
  int n_bits = 0;
  for (const auto& en : encodings) {
    n_bits += en.bit_precision;
  }
  std::vector<PhenConvert> conversions;
  for (const auto& dist : dists) {
    conversions.push_back(
        PhenConvert{encodings, Interval{lower(dist), upper(dist)}});
  }
  FirstPriceAuction auction(dists);
  std::vector<std::shared_ptr<multipop::AbstractSubGA<FirstPriceAuction>>>
      sub_gas;
  std::vector<std::shared_ptr<
      multipop::SinglePopulationGAToSubGAAdapter<FirstPriceAuction, Phen>>>
      ind_gas;

  int pop_size = 500;
  for (int i = 0; i < dists.size(); ++i) {
    auto ga = std::make_unique<SinglePopulationGA<Gen, Phen>>(
        BinaryGA<Phen>(conversions[i], pop_size, n_bits));

    auto subga = std::make_shared<
        multipop::SinglePopulationGAToSubGAAdapter<FirstPriceAuction, Phen>>(
        std::move(ga), i);
    sub_gas.push_back(subga);
    ind_gas.push_back(subga);
  }
  multipop::GA<FirstPriceAuction> ga(sub_gas, auction);
  ga.RunRound(500);
  std::vector<PhenotypeStrategy<Phen>> phens;
  for (int i = 0; i < dists.size(); ++i) {
    phens.push_back(ind_gas[i]->GetBestStrategy());
  }
  for (int i = std::floor(minValue); i <= std::ceil(maxValue); ++i) {
    std::cout << i << ": ";
    for (int player = 0; player < phens.size(); ++player) {
      if (i < std::floor(lower(dists[player])) ||
          i > std::ceil(upper(dists[player]))) {
        std::cout << "---" << '\t';
      } else {
        std::cout << phens[player].phenotype(i) << '\t';
      }
    }
    std::cout << '\n';
  }
  std::cout << "Expected Profits: ";
  for (auto strat : phens) {
    std::cout << "\t" << strat.fitness;
  }
  std::cout << std::endl;
}

void RunAllPay() {
  using PhenConvert = PhenConvertSort;
  std::vector<float> values{1, 1};
  std::vector<binary::Encoding> encodings(30, binary::Encoding{32, 0, 1});
  int n_bits = 0;
  for (const auto& en : encodings) {
    n_bits += en.bit_precision;
  }
  std::vector<PhenConvert> conversions;
  for (auto value : values) {
    conversions.push_back(PhenConvert{encodings, Interval{0, value}});
  }

  AllPayAuction auction(values);
  std::vector<std::shared_ptr<multipop::AbstractSubGA<AllPayAuction>>> sub_gas;
  std::vector<std::shared_ptr<
      multipop::SinglePopulationGAToSubGAAdapter<AllPayAuction, Phen>>>
      ind_gas;

  int pop_size = 100;
  for (int i = 0; i < values.size(); ++i) {
    auto ga = std::make_unique<SinglePopulationGA<Gen, Phen>>(
        BinaryGA<Phen>(conversions[i], pop_size, n_bits));
    auto subga = std::make_shared<
        multipop::SinglePopulationGAToSubGAAdapter<AllPayAuction, Phen>>(
            std::move(ga), i);
    sub_gas.push_back(subga);
    ind_gas.push_back(subga);
  }

  multipop::GA<AllPayAuction> ga(sub_gas, auction);
  ga.RunRound(10000);
  std::vector<PhenotypeStrategy<Phen>> phens;
  for (int i = 0; i < values.size(); ++i) {
    phens.push_back(ind_gas[i]->GetBestStrategy());
  }
  for (float i = 0; i <= 1; i += 0.01) {
    std::cout << i << ": ";
    for (int player = 0; player < phens.size(); ++player) {
      std::cout << phens[player].phenotype.GetBid(i) << '\t';
    }
    std::cout << '\n';
  }
  std::cout << "Expected Profits: ";
  for (auto strat : phens) {
    std::cout << "\t" << strat.fitness;
  }
  std::cout << std::endl;
}

int main(int argc, char** argv) {
  //RunAllPay();
  RunFirstPrice();
  return 0;
}
