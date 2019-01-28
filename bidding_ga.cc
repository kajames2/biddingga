#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/uniform.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "auctions/all_pay.h"
#include "auctions/common_value_signal.h"
#include "auctions/first_price.h"

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

#include "numericaldists/distribution.h"
#include "numericaldists/piecewise_linear.h"
#include "numericaldists/uneven_piecewise_linear.h"

using namespace genericga;
using namespace auctions;
using namespace numericaldists;
using namespace boost::math;

using Gen = binary::ByteArrayGenotype;
using Phen = PiecewiseLinear;

template <class Phen>
SinglePopulationGA<binary::ByteArrayGenotype, Phen> BinaryGA(
    std::function<Phen(const binary::ByteArrayGenotype&)> phen_conv,
    int pop_size, int n_bits,
    std::function<float(const Phen&)> fit = [](const Phen&) { return -1.0; }) {
  using Gen = binary::ByteArrayGenotype;
  int n_bytes = (n_bits + CHAR_BIT - 1) / CHAR_BIT;
  auto zero_genotype = Gen(std::vector<unsigned char>(n_bytes));
  std::vector<Gen> genes(pop_size, zero_genotype);
  Population<Gen, Phen> init_pop(phen_conv, fit, genes);

  auto children_fact = std::make_unique<ChildrenFactory<Gen>>(
      std::make_unique<binary::SinglePointCrossover>(),
      std::make_unique<binary::BitMutator>(3));

  return SinglePopulationGA<binary::ByteArrayGenotype, Phen>(
      std::move(init_pop), std::move(children_fact));
}

struct BinaryToPiecewiseRaw {
  PiecewiseLinear operator()(const Gen& gene) const {
    std::vector<float> ys = gene.ToFloatArray(nums_);
    return PiecewiseLinear(ys, interval_);
  }
  std::vector<binary::Encoding> nums_;
  Interval interval_;
};

struct BinaryToPiecewiseSort {
  PiecewiseLinear operator()(const Gen& gene) const {
    std::vector<float> ys = gene.ToFloatArray(nums_);
    std::sort(ys.begin(), ys.end());
    return PiecewiseLinear(ys, interval_);
  }
  std::vector<binary::Encoding> nums_;
  Interval interval_;
};

struct BinaryToPiecewiseFlooring {
  PiecewiseLinear operator()(const Gen& gene) const {
    std::vector<float> ys = gene.ToFloatArray(nums_);
    float floor = ys[0];
    for (auto& val : ys) {
      if (val < floor) {
        val = floor;
      }
      floor = val;
    }
    return PiecewiseLinear(ys, interval_);
  }
  std::vector<binary::Encoding> nums_;
  Interval interval_;
};

struct BinaryToPiecewiseCumulative {
  PiecewiseLinear operator()(const Gen& gene) const {
    std::vector<float> ys = gene.ToFloatArray(nums_);
    float tot = 0;
    for (auto& val : ys) {
      tot += val;
      if (tot > 1) {
        tot = 1;
      }
      val = tot;
    }
    return PiecewiseLinear(ys, interval_);
  }
  std::vector<binary::Encoding> nums_;
  Interval interval_;
};

struct BidFunctionGAConfiguration {
  Interval value_range = {0, 1};
  Interval bid_range = {0, 1};
  int n_strategies = 50;
  int n_children = 50;
  int n_segments = 30;
  int bit_precision = 32;
};

void RunFirstPrice() {
  using BinaryToPiecewise = BinaryToPiecewiseSort;

  std::vector<Distribution> dists{uniform_distribution<>(0, 100),
                                  uniform_distribution<>(0, 100),
                                  uniform_distribution<>(0, 100)};

  Interval value_range = {lower(dists), upper(dists)};
  std::vector<BidFunctionGAConfiguration> configs;
  for (const auto& dist : dists) {
    BidFunctionGAConfiguration config;
    config.value_range = {lower(dist), upper(dist)};
    config.bid_range = {0, value_range.max};
    config.n_strategies = 2000;
    config.n_children = 2000;
    config.n_segments = 30;
    config.bit_precision = 32;
    configs.push_back(config);
  }
  int n_rounds = 500;
  
  FirstPrice auction(dists);

  std::vector<std::shared_ptr<
      multipop::SinglePopulationGAToSubGAAdapter<FirstPrice, Phen>>>
      ind_gas;
  for (int i = 0; i < configs.size(); ++i) {
    const auto& config = configs[i];
    int n_bits = config.n_segments * config.bit_precision;
    binary::Encoding bid_enc{config.bit_precision, config.bid_range.min,
                             config.bid_range.max};
    std::vector<binary::Encoding> encodings(config.n_segments, bid_enc);
    BinaryToPiecewise conversion(
        BinaryToPiecewise{encodings, config.value_range});

    auto ga = std::make_unique<SinglePopulationGA<Gen, Phen>>(
        BinaryGA<Phen>(conversion, config.n_strategies, n_bits));

    auto subga = std::make_shared<
        multipop::SinglePopulationGAToSubGAAdapter<FirstPrice, Phen>>(
        std::move(ga), i);
    ind_gas.push_back(subga);
  }

  std::vector<std::shared_ptr<multipop::AbstractSubGA<FirstPrice>>> sub_gas(
      ind_gas.begin(), ind_gas.end());
  multipop::GA<FirstPrice> ga(sub_gas, auction);

  ga.RunRound(n_rounds);
  std::vector<PhenotypeStrategy<Phen>> phens;
  for (int i = 0; i < dists.size(); ++i) {
    phens.push_back(ind_gas[i]->GetBestStrategy());
  }
  for (int i = std::ceil(value_range.min); i <= std::floor(value_range.max);
       ++i) {
    std::cout << i << ": ";
    for (int player = 0; player < phens.size(); ++player) {
      if (!InInterval(configs[player].value_range, i)) {
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

void RunCommonValueSignal() {
  using BinaryToPiecewise = BinaryToPiecewiseSort;
  std::vector<int> n_draws{2, 2, 2, 2};
  float epsilon = 1;
  Interval discount_range{-1 * epsilon, 4 * epsilon};
  int n_lines = 30;
  int bit_precision = 32;
  std::vector<binary::Encoding> encodings(
      n_lines,
      binary::Encoding{bit_precision, discount_range.min, discount_range.max});
  int n_bits = 0;
  for (const auto& en : encodings) {
    n_bits += en.bit_precision;
  }

  std::vector<BinaryToPiecewise> conversions(
      n_draws.size(), BinaryToPiecewise{encodings, discount_range});

  CommonValueSignal auction(n_draws, epsilon, discount_range);
  std::vector<std::shared_ptr<multipop::AbstractSubGA<CommonValueSignal>>>
      sub_gas;
  std::vector<std::shared_ptr<
      multipop::SinglePopulationGAToSubGAAdapter<CommonValueSignal, Phen>>>
      ind_gas;

  int pop_size = 100;
  for (int i = 0; i < n_draws.size(); ++i) {
    auto ga = std::make_unique<SinglePopulationGA<Gen, Phen>>(
        BinaryGA<Phen>(conversions[i], pop_size, n_bits));

    auto subga = std::make_shared<
        multipop::SinglePopulationGAToSubGAAdapter<CommonValueSignal, Phen>>(
        std::move(ga), i);
    sub_gas.push_back(subga);
    ind_gas.push_back(subga);
  }
  multipop::GA<CommonValueSignal> ga(sub_gas, auction);
  for (int i = 0; i < 1; ++i) {
    ga.RunRound(2);
    std::cout << (i + 1) * 10 << "% complete" << std::endl;
  }

  std::vector<PhenotypeStrategy<Phen>> phens;
  for (int i = 0; i < n_draws.size(); ++i) {
    phens.push_back(ind_gas[i]->GetBestStrategy());
  }
  for (int i = 0; i <= 100; ++i) {
    float range = 2 * epsilon / 100 * i;
    std::cout << range << ": ";
    for (int player = 0; player < phens.size(); ++player) {
      std::cout << phens[player].phenotype(range) << '\t';
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
  using BinaryToPiecewise = BinaryToPiecewiseSort;
  std::vector<float> values{1, 1};
  std::vector<binary::Encoding> encodings(30, binary::Encoding{32, 0, 1});
  int n_bits = 0;
  for (const auto& en : encodings) {
    n_bits += en.bit_precision;
  }
  std::vector<BinaryToPiecewise> conversions;
  for (auto value : values) {
    conversions.push_back(BinaryToPiecewise{encodings, Interval{0, value}});
  }

  AllPay auction(values);
  std::vector<std::shared_ptr<multipop::AbstractSubGA<AllPay>>> sub_gas;
  std::vector<
      std::shared_ptr<multipop::SinglePopulationGAToSubGAAdapter<AllPay, Phen>>>
      ind_gas;

  int pop_size = 100;
  for (int i = 0; i < values.size(); ++i) {
    auto ga = std::make_unique<SinglePopulationGA<Gen, Phen>>(
        BinaryGA<Phen>(conversions[i], pop_size, n_bits));
    auto subga = std::make_shared<
        multipop::SinglePopulationGAToSubGAAdapter<AllPay, Phen>>(std::move(ga),
                                                                  i);
    sub_gas.push_back(subga);
    ind_gas.push_back(subga);
  }

  multipop::GA<AllPay> ga(sub_gas, auction);
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
  // RunAllPay();
  //RunFirstPrice();
  RunCommonValueSignal();
  return 0;
}
