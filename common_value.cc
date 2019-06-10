#include <variant>
#include <vector>

#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/uniform.hpp>

#include "auctions/all_pay.h"
#include "auctions/common_value_signal.h"
#include "auctions/common_value_signal_second.h"
#include "auctions/first_price.h"
#include "auctions/first_price_reverse.h"
#include "auctions/second_price.h"

#include "genericga/binary/bit_mutator.h"
#include "genericga/binary/byte_array_genotype.h"
#include "genericga/binary/encoding.h"
#include "genericga/binary/single_point_crossover.h"
#include "genericga/composite_ga.h"
#include "genericga/multipop/abstract_sub_ga.h"
#include "genericga/multipop/ga.h"
#include "genericga/multipop/sub_ga_adapter.h"
#include "genericga/selector/elitism_decorator.h"
#include "genericga/selector/keep_best.h"
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
  int id = 0;
  Interval value_range = {0, 1};
  Interval bid_range = {0, 1};
  int n_strategies = 50;
  int n_children = 50;
  int n_segments = 30;
  int bit_precision = 32;
};

struct BinaryToFloat {
  float operator()(const Gen& gene) const {
    return gene.ToFloatArray({num_})[0];
  }
  binary::FloatEncoding num_;
};

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

struct SeparateJoiner {
  PhenotypeStrategy<Phen> operator()(
      const std::vector<PhenotypeStrategy<Phen>>& strats) {
    int n_strats = strats.size();
    ArrayXd joined_xs(strat_size * n_strats);
    ArrayXd joined_ys(strat_size * n_strats);
    float fit = 0;
    for (int i = 0; i < n_strats; ++i) {
      fit += strats[i].fitness;
      joined_xs.segment(strat_size * i, strat_size) << strats[i].phenotype.xs;
      joined_ys.segment(strat_size * i, strat_size) << strats[i].phenotype.ys;
    }
    return PhenotypeStrategy<Phen>{{joined_xs, joined_ys}, fit};
  }
  int strat_size;
};

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
      std::make_unique<binary::BitMutator>(3),
      std::make_unique<selector::TournamentMixed>(2.2));

  return SinglePopulationGA<binary::ByteArrayGenotype, Phen>(
      std::move(init_pop), std::move(children_fact),
      std::make_unique<selector::ElitismDecorator>(
          std::make_unique<selector::TournamentMixed>(2.2), 1));
}

template <class Phen>
std::unique_ptr<SinglePopulationGA<Gen, Phen>> Make0DGA(
    BidFunctionGAConfiguration config) {
  int n_bits = config.bit_precision;
  binary::FloatEncoding bid_enc(static_cast<float>(config.bid_range.min),
                                static_cast<float>(config.bid_range.max),
                                config.bit_precision);
  binary::FloatEncoding encoding(bid_enc);
  std::function<Phen(const Gen&)> conversion = BinaryToFloat{encoding};
  return std::make_unique<SinglePopulationGA<Gen, Phen>>(
      BinaryGA<Phen>(conversion, config.n_strategies, n_bits));
}

template <class Phen>
std::unique_ptr<CompositeGA<Phen>> Make1DGA(BidFunctionGAConfiguration config,
                                            int n_composites) {
  assert(config.n_segments % n_composites == 0);
  int segs_per_comp = config.n_segments / n_composites;
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
std::shared_ptr<multipop::SubGAAdapter<Environment, Phen>> MakeSub1DGA(
    BidFunctionGAConfiguration config, int n_composites) {
  auto ga = Make1DGA<Phen>(std::move(config), n_composites);
  return std::make_shared<multipop::SubGAAdapter<Environment, Phen>>(
      std::move(ga), config.id);
}

template <class Environment, class Phen>
std::shared_ptr<multipop::SubGAAdapter<Environment, Phen>> MakeSub0DGA(
    BidFunctionGAConfiguration config) {
  std::shared_ptr<multipop::SubGAAdapter<Environment, Phen>> ind_ga;
  auto ga = Make0DGA<Phen>(std::move(config));
  auto subga = std::make_shared<multipop::SubGAAdapter<Environment, Phen>>(
      std::move(ga), config.id);
  return subga;
}

template <class Phen>
using subga_ptr =
    std::shared_ptr<multipop::SubGAAdapter<CommonValueSignalSecond, Phen>>;

int main(int argc, char** argv) {
  std::vector<int> n_draws;
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " Draws1, Draws2 ..." << std::endl;
    return 1;
  } else {
    for (int i = 1; i < argc; ++i) {
      n_draws.push_back(std::atoi(argv[i]));
    }
  }
  double epsilon = 500;
  Interval range_range{0, 2 * epsilon};
  Interval bid_range{-4 * epsilon, epsilon};

  int n_composites = 10;
  std::vector<BidFunctionGAConfiguration> configs;
  for (int i = 0; i < n_draws.size(); ++i) {
    BidFunctionGAConfiguration config;
    config.id = i;
    config.value_range = {0, 2 * epsilon};
    config.bid_range = {-4 * epsilon, epsilon};
    config.n_strategies = 500 / n_composites;
    config.n_children = 500 / n_composites;
    config.n_segments = 50;
    configs.push_back(std::move(config));
  }

  int n_rounds = 1000;
  CommonValueSignalSecond auction(n_draws, epsilon, {-4 * epsilon, epsilon});
  // FirstPrice auction(dists);

  std::vector<std::variant<subga_ptr<float>, subga_ptr<Phen>>> gas;
  for (int i = 0; i < n_draws.size(); ++i) {
    if (n_draws[i] == 1) {
      gas.push_back(
          MakeSub0DGA<CommonValueSignalSecond, float>(std::move(configs[i])));
    } else {
      gas.push_back(MakeSub1DGA<CommonValueSignalSecond, Phen>(std::move(configs[i]),
                                                         n_composites));
    }
  }

  std::vector<std::shared_ptr<multipop::AbstractSubGA<CommonValueSignalSecond>>>
      sub_gas;
  for (int i = 0; i < n_draws.size(); ++i) {
    if (n_draws[i] == 1) {
      sub_gas.push_back(std::get<subga_ptr<float>>(gas[i]));
    } else {
      sub_gas.push_back(std::get<subga_ptr<Phen>>(gas[i]));
    }
  }

  auto driver = multipop::GA<CommonValueSignalSecond>(sub_gas, auction);
  for (int n = 0; n < n_rounds; ++n) {
    driver.RunRound(1);
    for (int i = 0; i < gas.size(); ++i) {
      if (n_draws[i] == 1) {
        std::cout
            << std::get<subga_ptr<float>>(gas[i])->GetBestStrategy().phenotype
            << std::endl;
        auction.AcceptStrategy(
            std::get<subga_ptr<float>>(gas[i])->GetBestStrategy().phenotype, i);
      } else {
        std::cout
            << std::get<subga_ptr<Phen>>(gas[i])->GetBestStrategy().phenotype
            << std::endl;
        auction.AcceptStrategy(
            std::get<subga_ptr<Phen>>(gas[i])->GetBestStrategy().phenotype, i);
      }
    }
    for (int i = 0; i < gas.size(); ++i) {
      if (n_draws[i] == 1) {
        std::cout << auction.GetFitness(
            std::get<subga_ptr<float>>(gas[i])->GetBestStrategy().phenotype, i);
        if (i != gas.size() - 1) {
          std::cout << ",";
        }
      } else {
        std::cout << auction.GetFitness(
            std::get<subga_ptr<Phen>>(gas[i])->GetBestStrategy().phenotype, i);
        if (i != gas.size() - 1) {
          std::cout << ",";
        }
      }
    }
    std::cout << std::endl;
  }
}
