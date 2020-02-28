#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/uniform.hpp>

#include "auctions/common_value_signal_endpoints.h"

#include "biddingga/configuration.h"
#include "biddingga/initializers_1d.h"
#include "biddingga/initializers_2d.h"
#include "genericga/binary/bit_mutator.h"
#include "genericga/binary/byte_array_genotype.h"
#include "genericga/binary/encoding.h"
#include "genericga/binary/single_point_crossover.h"
#include "genericga/composite_ga.h"
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
#include "numericaldists/grid.h"
#include "numericaldists/scatter.h"

#include <cstdlib>
#include <eigen3/Eigen/Core>
#include <iostream>
#include <variant>
#include <vector>

using namespace genericga;
using namespace auctions;
using namespace numericaldists;
using namespace boost::math;
using namespace Eigen;
using namespace biddingga;

template <class Environment>
std::shared_ptr<multipop::SubGAAdapter<Environment, Scatter>> MakeSub1DGA(
    Configuration1D config) {
  int id = config.id;
  auto ga = MakeGAComposite(std::move(config));
  return std::make_shared<multipop::SubGAAdapter<Environment, Scatter>>(
      std::move(ga), id);
}

template <class Environment>
std::shared_ptr<multipop::SubGAAdapter<Environment, Grid>> MakeSub2DGA(
    Configuration2D config) {
  int id = config.id;
  auto ga = MakeGAComposite(std::move(config));
  return std::make_shared<multipop::SubGAAdapter<Environment, Grid>>(
      std::move(ga), id);
}
template <class Phen>
using subga_ptr =
    std::shared_ptr<multipop::SubGAAdapter<CommonValueSignalEndpoints, Phen>>;

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
  double min_value = 500;
  double max_value = 9500;
  Interval midpoint_range{min_value - epsilon, max_value + epsilon};
  Interval uncertainty_range{0, 2 * epsilon};
  uniform_distribution<> value_dist =
      uniform_distribution<>(min_value, max_value);
  uniform_distribution<> error_dist = uniform_distribution<>(-epsilon, epsilon);
  CommonValueSignalEndpoints auction(value_dist, error_dist, n_draws);

  std::vector<std::variant<subga_ptr<Scatter>, subga_ptr<Grid>>> gas;
  std::vector<
      std::shared_ptr<multipop::AbstractSubGA<CommonValueSignalEndpoints>>>
      sub_gas;
  for (int i = 0; i < n_draws.size(); ++i) {
    if (n_draws[i] == 1) {
      Configuration1D config;
      config.id = i;
      config.x_range = midpoint_range;
      config.y_range = midpoint_range;
      config.n_strategies = 100;
      config.n_children = 100;
      gas.push_back(MakeSub1DGA<CommonValueSignalEndpoints>(std::move(config)));
      sub_gas.push_back(std::get<subga_ptr<Scatter>>(gas[i]));
    } else {
      Configuration2D config;
      config.id = i;
      config.x_range = midpoint_range;
      config.y_range = uncertainty_range;
      config.z_range = midpoint_range;
      gas.push_back(MakeSub2DGA<CommonValueSignalEndpoints>(std::move(config)));
      sub_gas.push_back(std::get<subga_ptr<Grid>>(gas[i]));
    }
  }
  int n_rounds = 3000;
  int output_frequency = 30;
  auto driver = multipop::GA<CommonValueSignalEndpoints>(sub_gas, auction);
  for (int n = 0; n < n_rounds; ++n) {
    driver.RunRound(1);
    if (n % output_frequency == 0) {
      for (int i = 0; i < gas.size(); ++i) {
        if (n_draws[i] == 1) {
          std::cout << std::get<subga_ptr<Scatter>>(gas[i])
                           ->GetBestStrategy()
                           .phenotype
                    << std::endl;
          auction.AcceptStrategy(
              std::get<subga_ptr<Scatter>>(gas[i])->GetBestStrategy().phenotype,
              i);
        } else {
          std::cout
              << std::get<subga_ptr<Grid>>(gas[i])->GetBestStrategy().phenotype
              << std::endl;
          auction.AcceptStrategy(
              std::get<subga_ptr<Grid>>(gas[i])->GetBestStrategy().phenotype,
              i);
        }
      }
      for (int i = 0; i < gas.size(); ++i) {
        if (n_draws[i] == 1) {
          std::cout << auction.GetFitness(
              std::get<subga_ptr<Scatter>>(gas[i])->GetBestStrategy().phenotype,
              i);
          if (i != gas.size() - 1) {
            std::cout << ",";
          }
        } else {
          std::cout << auction.GetFitness(
              std::get<subga_ptr<Grid>>(gas[i])->GetBestStrategy().phenotype,
              i);
          if (i != gas.size() - 1) {
            std::cout << ",";
          }
        }
      }
      std::cout << std::endl;
    }
  }
  return 0;
}
