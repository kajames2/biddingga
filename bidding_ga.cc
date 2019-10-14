#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/uniform.hpp>

#include "auctions/all_pay.h"
#include "auctions/first_price.h"
#include "auctions/first_price_reverse.h"
#include "auctions/second_price.h"
#include "auctions/common_value_endpoints.h"
#include "auctions/common_value_endpoints2.h"

#include "genericga/binary/bit_mutator.h"
#include "genericga/binary/byte_array_genotype.h"
#include "genericga/binary/encoding.h"
#include "genericga/binary/single_point_crossover.h"
#include "genericga/composite_ga.h"
#include "genericga/multipop/ga.h"
#include "genericga/multipop/sub_ga_adapter.h"
#include "genericga/selector/elitism_decorator.h"
#include "genericga/selector/keep_best.h"
#include "genericga/selector/keep_commonest.h"
#include "genericga/selector/ranked_weighted.h"
#include "genericga/selector/tournament.h"
#include "genericga/selector/tournament_mixed.h"
#include "genericga/selector/tournament_poisson.h"
#include "genericga/single_population_ga.h"

#include "numericaldists/distribution.h"

#include <eigen3/Eigen/Core>

using namespace genericga;
using namespace auctions;
using namespace numericaldists;
using namespace boost::math;
using namespace Eigen;

using Gen = binary::ByteArrayGenotype;
using Phen = Scatter;

struct BidFunctionGAConfiguration {
  int id = -1;
  Interval value_range = {0, 1};
  Interval bid_range = {0, 1};
  int n_composites = 5;
  int n_strategies = 100;
  int n_children = 100;
  int n_segments = 100;
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
      std::make_unique<binary::BitMutator>(2),
      std::make_unique<selector::TournamentMixed>(1));

  return SinglePopulationGA<binary::ByteArrayGenotype, Phen>(
      std::move(init_pop), std::move(children_fact),
      std::make_unique<selector::ElitismDecorator>(
          std::make_unique<selector::TournamentMixed>(2.2), 2));
}

struct SortDecorator {
  Scatter operator()(const Gen& gene) {
    Scatter scatter = phen_conv(gene);
    std::sort(scatter.ys.data(), scatter.ys.data() + scatter.ys.size());
    return scatter;
  }
  std::function<Phen(const Gen& gene)> phen_conv;
};

struct ReverseSortDecorator {
  Scatter operator()(const Gen& gene) {
    Scatter scatter = phen_conv(gene);
    std::sort(scatter.ys.data(), scatter.ys.data() + scatter.ys.size());
    std::reverse(scatter.ys.data(), scatter.ys.data() + scatter.ys.size());
    return scatter;
  }
  std::function<Phen(const Gen&)> phen_conv;
};

struct BinaryToScatter {
  Scatter operator()(const Gen& gene) const {
    ArrayXf ys = gene.ToEigenFloatArray(nums_);
    return {vals_, ys.cast<double>()};
  }
  ArrayXd vals_;
  std::vector<binary::FloatEncoding> nums_;
};

struct JoinerSortDecorator {
  PhenotypeStrategy<Phen> operator()(
      const std::vector<PhenotypeStrategy<Phen>>& strats) {
    auto strat = joiner(strats);
    std::sort(strat.phenotype.ys.data(),
              strat.phenotype.ys.data() + strat.phenotype.ys.size());
    return strat;
  }
  std::function<PhenotypeStrategy<Phen>(
      const std::vector<PhenotypeStrategy<Phen>>&)>
      joiner;
};

struct JoinerReverseSortDecorator {
  PhenotypeStrategy<Phen> operator()(
      const std::vector<PhenotypeStrategy<Phen>>& strats) {
    auto strat = joiner(strats);
    std::sort(strat.phenotype.ys.data(),
              strat.phenotype.ys.data() + strat.phenotype.ys.size());
    std::reverse(strat.phenotype.ys.data(),
                 strat.phenotype.ys.data() + strat.phenotype.ys.size());
    return strat;
  }
  std::function<PhenotypeStrategy<Phen>(
      const std::vector<PhenotypeStrategy<Phen>>&)>
      joiner;
};

// struct SeparateJoiner {
//   PhenotypeStrategy<Phen> operator()(
//       const std::vector<PhenotypeStrategy<Phen>>& strats) {
//     int n_strats = strats.size();
//     ArrayXd joined_xs(strat_size * n_strats);
//     ArrayXd joined_ys(strat_size * n_strats);
//     float fit = 0;
//     for (int i = 0; i < n_strats; ++i) {
//       fit += strats[i].fitness;
//       joined_xs.segment(strat_size * i, strat_size) <<
//       strats[i].phenotype.xs; joined_ys.segment(strat_size * i, strat_size)
//       << strats[i].phenotype.ys;
//     }
//     return PhenotypeStrategy<Phen>{{joined_xs, joined_ys}, fit};
//   }
//   int strat_size;
// };

struct MedianJoiner {
  PhenotypeStrategy<Phen> operator()(
      const std::vector<PhenotypeStrategy<Phen>>& strats) {
    int n_strats = strats.size();
    ArrayXd joined_xs(n_strats * (strat_size - 1) + 1);
    ArrayXd joined_ys(n_strats * (strat_size - 1) + 1);
    float fit = 0;
    joined_xs(0) = strats[0].phenotype.xs(0);
    joined_ys(0) = strats[0].phenotype.ys(0);
    for (int i = 0; i < n_strats - 1; ++i) {
      fit += strats[i].fitness;
      joined_xs.segment(i * (strat_size - 1) + 1, strat_size - 1)
          << strats[i].phenotype.xs.segment(1, strat_size - 1);
      joined_ys.segment(i * (strat_size - 1) + 1, strat_size - 2)
          << strats[i].phenotype.ys.segment(1, strat_size - 2);
      joined_ys((i + 1) * (strat_size - 1)) =
          (strats[i].phenotype.ys(strat_size - 1) +
           strats[i + 1].phenotype.ys(0)) /
          2;
    }
    fit += strats[n_strats - 1].fitness;
    joined_xs.segment((n_strats - 1) * (strat_size - 1) + 1, strat_size - 1)
        << strats[strats.size() - 1].phenotype.xs.segment(1, strat_size - 1);
    joined_ys.segment((n_strats - 1) * (strat_size - 1) + 1, strat_size - 1)
        << strats[strats.size() - 1].phenotype.ys.segment(1, strat_size - 1);
    return PhenotypeStrategy<Phen>{{joined_xs, joined_ys}, fit};
  }
  int strat_size;
};

template <class Phen>
std::unique_ptr<CompositeGA<Phen>> MakeGAComposite(
    BidFunctionGAConfiguration config) {
  assert(config.n_segments % config.n_composites == 0);
  int segs_per_comp = config.n_segments / config.n_composites;
  int n_floats_per_comp = segs_per_comp + 1;
  int n_bits = n_floats_per_comp * config.bit_precision;
  binary::FloatEncoding bid_enc(static_cast<float>(config.bid_range.min),
                                static_cast<float>(config.bid_range.max),
                                config.bit_precision);
  std::vector<binary::FloatEncoding> encodings(n_floats_per_comp, bid_enc);
  ArrayXd vals = ArrayXd::LinSpaced(
      config.n_segments + 1, config.value_range.min, config.value_range.max);
  std::vector<std::shared_ptr<AbstractSinglePopulationGA<Phen>>> gas;
  for (int i = 0; i < config.n_segments; i += segs_per_comp) {
    ArrayXd sub_vals = vals.segment(i, n_floats_per_comp);
    std::function<Phen(const Gen&)> conversion =
        BinaryToScatter{sub_vals, encodings};
    gas.push_back(std::make_shared<SinglePopulationGA<Gen, Phen>>(
        BinaryGA<Phen>(conversion, config.n_strategies, n_bits)));
  }
  return std::make_unique<CompositeGA<Phen>>(
      gas, MedianJoiner{n_floats_per_comp});
}

template <class Environment, class Phen>
std::vector<std::shared_ptr<multipop::SubGAAdapter<Environment, Phen>>>
MakeSubGAs(std::vector<BidFunctionGAConfiguration> configs) {
  std::vector<std::shared_ptr<multipop::SubGAAdapter<Environment, Phen>>>
      ind_gas;
  for (int i = 0; i < configs.size(); ++i) {
    auto ga = MakeGAComposite<Phen>(std::move(configs[i]));
    auto subga = std::make_shared<multipop::SubGAAdapter<Environment, Phen>>(
        std::move(ga), i, 0, std::make_unique<selector::KeepBest>());
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

template <class Environment, class Phen>
void RunAndOutput(
    multipop::GA<Environment>& driver,
    std::vector<std::shared_ptr<multipop::SubGAAdapter<Environment, Phen>>>&
        gas,
    Environment& env, int n_rounds, int output_frequency) {
  for (int n = 0; n < n_rounds / output_frequency; ++n) {
    driver.RunRound(output_frequency);
    for (int j = 0; j < gas.size(); ++j) {
      ArrayXd xvals = gas[j]->GetBestStrategy().phenotype.xs.transpose();
      for (int i = 0; i < xvals.size(); ++i) {
        std::cout << xvals(i);
        if (i != xvals.size() - 1) {
          std::cout << ",";
        }
      }
      std::cout << std::endl;
      ArrayXd vals = gas[j]->GetBestStrategy().phenotype.ys.transpose();
      for (int i = 0; i < vals.size(); ++i) {
        std::cout << vals(i);
        if (i != vals.size() - 1) {
          std::cout << ",";
        }
      }
      std::cout << std::endl;
      env.AcceptStrategy(gas[j]->GetBestStrategy().phenotype, j);
    }
    for (int j = 0; j < gas.size(); ++j) {
      std::cout << env.GetFitness(gas[j]->GetBestStrategy().phenotype, j);
      if (j != gas.size() - 1) {
        std::cout << ",";
      }
    }
    std::cout << std::endl;
    // for (int j = 0; j < gas.size(); ++j) {
    //   std::cout << env.GetRevenue(gas[j]->GetBestStrategy().phenotype, j);
    //   if (j != gas.size() - 1) {
    //     std::cout << ",";
    //   }
    // }
    // std::cout << std::endl;
    // for (int j = 0; j < gas.size(); ++j) {
    //   std::cout << env.GetValue(gas[j]->GetBestStrategy().phenotype, j);
    //   if (j != gas.size() - 1) {
    //     std::cout << ",";
    //   }
    // }
    // std::cout << std::endl;
  }
}

int main(int argc, char** argv) {
  // std::vector<Distribution> dists{
  //       uniform_distribution<>(0, 1),
  //       uniform_distribution<>(0, 1)};
  // std::vector<Distribution> dists{uniform_distribution<>(0, 1),
  //                                uniform_distribution<>(0, 1),
  //                                normal_distribution<>(0.7, 0.2)};

  // Interval value_range = {lower(dists), upper(dists)};
  // std::vector<BidFunctionGAConfiguration> configs;
  // for (int i = 0; i < dists.size(); ++i) {
  //   BidFunctionGAConfiguration config;
  //   config.value_range = {lower(dists[i]), upper(dists[i])};
  //   config.bid_range = {0, value_range.max};
  //   configs.push_back(config);
  // }

  int n_players = 4;
  Distribution value_dist = uniform_distribution<>(500, 9500);
  Distribution error_dist = uniform_distribution<>(-500, 500);
  std::vector<BidFunctionGAConfiguration> configs;
  for (int i = 0; i < n_players; ++i) {
    BidFunctionGAConfiguration config;
    config.value_range = {lower(value_dist) + lower(error_dist), upper(value_dist) + upper(error_dist)};
    config.bid_range = {lower(value_dist) + lower(error_dist), upper(value_dist) + upper(error_dist)};
    configs.push_back(config);
  }

  
  using Auction = CommonValueEndpoints;
  Auction auction(n_players, value_dist, error_dist);
  //Auction auction(n_players, value_dist, error_dist);
  auto gas = MakeSubGAs<Auction, Scatter>(configs);
  auto driver = MakeMultipopDriver<Auction, Scatter>(gas, auction);

  int n_rounds = 1000;
  int output_frequency = 10;
  RunAndOutput(driver, gas, auction, n_rounds, output_frequency);

  return 0;
}
