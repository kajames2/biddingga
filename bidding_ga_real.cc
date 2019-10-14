#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/uniform.hpp>

#include "auctions/all_pay.h"
#include "auctions/first_price.h"
#include "auctions/second_price.h"

#include "genericga/composite_ga.h"
#include "genericga/multipop/ga.h"
#include "genericga/multipop/sub_ga_adapter.h"
#include "genericga/real/modify_mutation.h"
#include "genericga/real/single_point_crossover.h"
#include "genericga/selector/elitism_decorator.h"
#include "genericga/selector/keep_best.h"
#include "genericga/selector/ranked_weighted.h"
#include "genericga/selector/tournament.h"
#include "genericga/single_population_ga.h"

#include "numericaldists/distribution.h"

#include <eigen3/Eigen/Core>

using namespace genericga;
using namespace auctions;
using namespace numericaldists;
using namespace boost::math;
using namespace Eigen;

using Gen = std::vector<double>;
using Phen = Scatter;

struct BidFunctionGAConfiguration {
  Interval value_range = {0, 1};
  Interval bid_range = {0, 1};
  int n_strategies = 50;
  int n_children = 50;
  int n_segments = 30;
};

template <class Phen>
SinglePopulationGA<Gen, Phen> GA1DFunc(
    std::function<Phen(const Gen&)> phen_conv, int pop_size, int n_segments,
    double min, double max,
    std::function<std::vector<float>(const std::vector<Phen>&)> fit =
        [](const std::vector<Phen>& phens) {
          return std::vector<float>(phens.size(), -1.0);
        }) {
  std::vector<Gen> genes;
  auto generator = std::mt19937(std::random_device()());
  std::uniform_real_distribution<> dist(min, max);
  for (int i = 0; i < pop_size; ++i) {
    std::vector<double> rand_gene(n_segments + 1);
    for (double& r : rand_gene) {
      r = dist(generator);
    }
    genes.emplace_back(rand_gene);
  }
  Population<Gen, Phen> init_pop(phen_conv, fit, genes);
  auto children_fact = std::make_unique<ChildrenFactory<Gen>>(
      std::make_unique<real::SinglePointCrossover>(),
      std::make_unique<real::ModifyMutation>(3, (max - min) / 5., min, max),
      std::make_unique<selector::RankedWeighted>(0.6));

  return SinglePopulationGA<Gen, Phen>(
      std::move(init_pop), std::move(children_fact),
      std::make_unique<selector::ElitismDecorator>(
          std::make_unique<selector::RankedWeighted>(0.6), 1));
}

struct VecToScatterSort {
  Scatter operator()(Gen gene) const {
    Eigen::ArrayXd ys(gene.size());
    for (int i = 0; i < ys.size(); ++i) {
      ys(i) = gene[i];
    }
    std::sort(ys.data(), ys.data() + ys.size());
    return {vals_, ys};
  }
  ArrayXd vals_;
};

template <class VecToPiecewise, class Phen>
std::unique_ptr<SinglePopulationGA<Gen, Phen>> MakeGA(
    BidFunctionGAConfiguration config) {
  int gene_size = config.n_segments + 1;
  ArrayXd vals = ArrayXd::LinSpaced(
      config.n_segments + 1, config.value_range.min, config.value_range.max);
  VecToPiecewise conversion(VecToPiecewise{vals});

  return std::make_unique<SinglePopulationGA<Gen, Phen>>(
      GA1DFunc<Phen>(conversion, config.n_strategies, config.n_segments,
                     config.bid_range.min, config.bid_range.max));
}

PhenotypeStrategy<Phen> SeparateJoiner(
    const std::vector<PhenotypeStrategy<Phen>>& strats) {
  ArrayXd joined_xs(2 * strats.size());
  ArrayXd joined_ys(2 * strats.size());
  float fit = 0;
  for (int i = 0; i < strats.size(); ++i) {
    fit += strats[i].fitness;
    joined_xs.segment(2 * i, 2) << strats[i].phenotype.xs;
    joined_ys.segment(2 * i, 2) << strats[i].phenotype.ys;
  }
  return PhenotypeStrategy<Phen>{{joined_xs, joined_ys}, fit};
}

PhenotypeStrategy<Phen> MedianJoiner(
    const std::vector<PhenotypeStrategy<Phen>>& strats) {
  ArrayXd joined_xs(strats.size() + 1);
  ArrayXd joined_ys(strats.size() + 1);
  float fit = 0;
  joined_xs.segment(0, 1) << strats[0].phenotype.xs(0);
  joined_ys.segment(0, 1) << strats[0].phenotype.ys(0);
  for (int i = 1; i < strats.size(); ++i) {
    fit += strats[i].fitness;
    joined_xs.segment(i, 1) << strats[i].phenotype.xs(0);
    joined_ys.segment(i, 1)
        << (strats[i - 1].phenotype.ys(1) + strats[i].phenotype.ys(0)) / 2;
  }
  joined_xs.segment(strats.size(), 1)
      << strats[strats.size() - 1].phenotype.xs(1);
  joined_ys.segment(strats.size(), 1)
      << strats[strats.size() - 1].phenotype.ys(1);
  return PhenotypeStrategy<Phen>{{joined_xs, joined_ys}, fit};
}

template <class VecToPiecewise, class Phen>
std::unique_ptr<CompositeGA<Phen>> MakeGAComposite(
    BidFunctionGAConfiguration config) {
  int gene_size = 2;
  ArrayXd vals = ArrayXd::LinSpaced(
      config.n_segments + 1, config.value_range.min, config.value_range.max);
  std::vector<std::shared_ptr<AbstractSinglePopulationGA<Phen>>> gas;
  for (int i = 0; i < config.n_segments; ++i) {
    ArrayXd sub_vals = vals.segment(i, 2);
    VecToPiecewise conversion(VecToPiecewise{sub_vals});
    gas.push_back(std::make_shared<SinglePopulationGA<Gen, Phen>>(
        GA1DFunc<Phen>(conversion, config.n_strategies, 1, config.bid_range.min,
                       config.bid_range.max)));
  }
  return std::make_unique<CompositeGA<Phen>>(gas, MedianJoiner);
}

template <class Environment, class Phen,
          class VecToPiecewise = VecToScatterSort>
std::vector<std::shared_ptr<multipop::SubGAAdapter<Environment, Phen>>>
MakeSubGAs(std::vector<BidFunctionGAConfiguration> configs) {
  std::vector<std::shared_ptr<multipop::SubGAAdapter<Environment, Phen>>>
      ind_gas;
  for (int i = 0; i < configs.size(); ++i) {
    auto ga = MakeGAComposite<VecToPiecewise, Phen>(std::move(configs[i]));
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
  return multipop::GA<Environment>(sub_gas, env, 10, 200);
}

void RunWinnerPay() {
  std::vector<Distribution> dists{uniform_distribution<>(0, 1),
                                  uniform_distribution<>(0, 1),
                                  uniform_distribution<>(0, 1)};

  Interval value_range = {lower(dists), upper(dists)};
  std::vector<BidFunctionGAConfiguration> configs;
  for (int i = 0; i < dists.size(); ++i) {
    BidFunctionGAConfiguration config;
    config.value_range = {lower(dists[i]), upper(dists[i])};
    config.bid_range = {0, value_range.max};
    config.n_strategies = 20;
    config.n_children = 20;
    config.n_segments = 100;
    configs.push_back(std::move(config));
  }
  int n_rounds = 1000;
  AllPay auction(dists);
  auto gas = MakeSubGAs<AllPay, Phen>(std::move(configs));
  auto driver = MakeMultipopDriver<AllPay, Phen>(gas, auction);
  ArrayXd vals = gas[0]->GetBestStrategy().phenotype.xs.transpose();
  for (int n = 0; n < n_rounds; ++n) {
    driver.RunRound(1);
    for (int j = 0; j < gas.size(); ++j) {
      ArrayXd xvals = gas[j]->GetBestStrategy().phenotype.xs.transpose();
      for (int i = 0; i < xvals.size(); ++i) {
        std::cout << xvals(i) << ",";
      }
      std::cout << std::endl;
      ArrayXd vals = gas[j]->GetBestStrategy().phenotype.ys.transpose();
      for (int i = 0; i < vals.size(); ++i) {
        std::cout << vals(i) << ",";
      }
      std::cout << std::endl;
    }
  }
}

int main(int argc, char** argv) {
  RunWinnerPay();
  return 0;
}
