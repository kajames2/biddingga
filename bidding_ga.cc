#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/uniform.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "auctions/all_pay.h"
#include "auctions/common_value_signal.h"
#include "auctions/winner_pay.h"

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

struct BidFunctionGAConfiguration {
  Interval value_range = {0, 1};
  Interval bid_range = {0, 1};
  int n_strategies = 50;
  int n_children = 50;
  int n_segments = 30;
  int bit_precision = 32;
};

template <class Phen>
SinglePopulationGA<binary::ByteArrayGenotype, Phen> BinaryGA(
    std::function<Phen(const binary::ByteArrayGenotype&)> phen_conv,
    int pop_size, int n_bits,
    std::function<std::vector<float>(const std::vector<Phen>&)> fit =
        [](const std::vector<Phen>& phens) {
          return std::vector<float>(phens.size(), -1.0);
        }) {
  using Gen = binary::ByteArrayGenotype;
  int n_bytes = (n_bits + CHAR_BIT - 1) / CHAR_BIT;
  std::vector<Gen> genes;
  auto generator = std::mt19937(std::random_device()());
  std::uniform_int_distribution<int> dist(0, UCHAR_MAX);
  for (int i = 0; i < pop_size; ++i) {
    std::vector<unsigned char> rand_gene(n_bytes);
    for (int j = 0; j < n_bytes; ++j) {
      rand_gene[j] = static_cast<unsigned char>(dist(generator));
    }
    genes.emplace_back(rand_gene);
  }
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

template <class BinaryToPiecewise, class Phen>
std::unique_ptr<SinglePopulationGA<binary::ByteArrayGenotype, Phen>> MakeGA(
    BidFunctionGAConfiguration config) {
  int n_bits = config.n_segments * config.bit_precision;
  binary::Encoding bid_enc{config.bit_precision, config.bid_range.min,
                           config.bid_range.max};
  std::vector<binary::Encoding> encodings(config.n_segments, bid_enc);
  BinaryToPiecewise conversion(
      BinaryToPiecewise{encodings, config.value_range});

  return std::make_unique<SinglePopulationGA<Gen, Phen>>(
      BinaryGA<Phen>(conversion, config.n_strategies, n_bits));
}

template <class Environment, class Phen,
          class BinaryToPiecewise = BinaryToPiecewiseSort>
std::vector<std::shared_ptr<multipop::SubGAAdapter<Environment, Phen>>>
MakeSubGAs(std::vector<BidFunctionGAConfiguration> configs) {
  std::vector<std::shared_ptr<multipop::SubGAAdapter<Environment, Phen>>>
      ind_gas;
  for (int i = 0; i < configs.size(); ++i) {
    auto ga = MakeGA<BinaryToPiecewise, Phen>(std::move(configs[i]));
    auto subga = std::make_shared<multipop::SubGAAdapter<Environment, Phen>>(
        std::move(ga), i);
    ind_gas.push_back(subga);
  }
  return ind_gas;
}

template <class Environment, class Phen>
multipop::GA<Environment> MakeMultipopDriver(
    const std::vector<
        std::shared_ptr<multipop::SubGAAdapter<Environment, Phen>>>& gas,
    Environment env) {
  std::vector<std::shared_ptr<multipop::AbstractSubGA<Environment>>> sub_gas(
      gas.begin(), gas.end());
  return multipop::GA<Environment>(sub_gas, env);
}

void RunWinnerPay() {
  std::vector<Distribution> dists{uniform_distribution<>(0, 1),
                                  uniform_distribution<>(0, 1),
                                  uniform_distribution<>(0, 1),
                                  uniform_distribution<>(0, 1),
                                  uniform_distribution<>(0, 1)};

  Interval value_range = {lower(dists), upper(dists)};
  std::vector<BidFunctionGAConfiguration> configs;
  // for (const auto& dist : dists) {
  for (int i = 0; i < 1; ++i) {
    BidFunctionGAConfiguration config;
    config.value_range = {lower(dists[i]), upper(dists[i])};
    config.bid_range = {0, 10 * value_range.max};
    config.n_strategies = 2000;
    config.n_children = 2000;
    config.n_segments = 100;
    config.bit_precision = 32;
    configs.push_back(std::move(config));
  }
  int n_rounds = 1000;
  WinnerPay auction(dists, 3, 10);
  for (int i = 0; i < dists.size(); ++i) {
    auction.AcceptStrategy(numericaldists::PiecewiseLinear{{0, 4/3.}, {0, 1}}, i);
  }
  std::cout << auction.GetFitness(numericaldists::PiecewiseLinear{{0, 4/3.}, {0, 1}}, 0) << std::endl;
  auto gas = MakeSubGAs<WinnerPay, Phen>(std::move(configs));
  auto driver = MakeMultipopDriver<WinnerPay, Phen>(gas, auction);
  driver.RunRound(n_rounds);

  std::vector<PhenotypeStrategy<Phen>> phens;
  for (int i = 0; i < gas.size(); ++i) {
    phens.push_back(gas[i]->GetBestStrategy());
  }
  
  for (int i = 0; i <= 100; ++i) {
    float value = GetSpan(value_range) / 100 * i + value_range.min;
    std::cout << std::setw(4) << std::setprecision(2) << value << ",";
    for (int player = 0; player < phens.size(); ++player) {
      if (!InInterval({lower(dists[player]), upper(dists[player])}, value)) {
        std::cout << "---";
      } else {
        std::cout << std::setw(10) << std::setprecision(3) << phens[player].phenotype(value);
      }
      if (player < phens.size() -1) {
        std::cout << ", ";
      }
    }
    std::cout << '\n';
  }
  std::cout << "Expected Profits, ";
  for (auto strat : phens) {
    std::cout << "\t" << strat.fitness;
  }
  std::cout << std::endl;
}

void RunCommonValueSignal() {
  std::vector<int> n_draws{2, 2, 2, 2};
  float epsilon = 1;
  Interval range_range{0, 2 * epsilon};
  Interval discount_range{-1 * epsilon, 4 * epsilon};

  std::vector<BidFunctionGAConfiguration> configs;
  for (int i = 0; i < n_draws.size(); ++i) {
    BidFunctionGAConfiguration config;
    config.value_range = range_range;
    config.bid_range = discount_range;
    config.n_strategies = 1000;
    config.n_children = 1000;
    config.n_segments = 30;
    config.bit_precision = 32;
    configs.push_back(std::move(config));
  }

  CommonValueSignal auction(n_draws, epsilon, discount_range);
  auto gas = MakeSubGAs<CommonValueSignal, Phen>(std::move(configs));
  auto driver = MakeMultipopDriver<CommonValueSignal, Phen>(gas, auction);

  for (int i = 0; i < 10; ++i) {
    driver.RunRound(10);
    if (i > 0) {
      std::cout << '\r';
    }
    std::cout << (i + 1) * 10 << "% complete" << std::flush;
  }
  std::cout << std::endl;
  std::vector<PhenotypeStrategy<Phen>> phens;
  for (int i = 0; i < n_draws.size(); ++i) {
    phens.push_back(gas[i]->GetBestStrategy());
  }

  std::cout << std::fixed;
  for (int i = 0; i <= 100; ++i) {
    float value = GetSpan(range_range) / 100 * i + range_range.min;
    std::cout << std::setprecision(2) << std::setw(4) << value << ", ";
    for (int player = 0; player < phens.size(); ++player) {
      std::cout << std::setprecision(3) << std::setw(5)
                << phens[player].phenotype(value);
      if (player < phens.size() - 1) {
        std::cout << ", ";
      }
    }
    std::cout << '\n';
  }
  std::cout << "Expected Profits";
  for (auto strat : phens) {
    std::cout << ", " << std::setw(6) << std::setprecision(3) << strat.fitness;
  }
  std::cout << std::endl;
}

void RunAllPay() {
  std::vector<float> values{1, 1};
  std::vector<BidFunctionGAConfiguration> configs;
  for (float value : values) {
    BidFunctionGAConfiguration config;
    config.value_range = {0, value};
    config.bid_range = {0, 1};
    config.n_strategies = 2000;
    config.n_children = 2000;
    config.n_segments = 30;
    config.bit_precision = 32;
    configs.push_back(std::move(config));
  }
  int n_rounds = 10000;
  AllPay auction(values);

  auto gas = MakeSubGAs<AllPay, Phen>(std::move(configs));
  auto driver = MakeMultipopDriver<AllPay, Phen>(gas, auction);
  driver.RunRound(n_rounds);

  std::vector<PhenotypeStrategy<Phen>> phens;
  for (int i = 0; i < values.size(); ++i) {
    phens.push_back(gas[i]->GetBestStrategy());
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
  RunAllPay();
  // RunWinnerPay();
  // RunCommonValueSignal();
  return 0;
}
